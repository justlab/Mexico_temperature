library(data.table)
library(FNN)
library(lme4)
library(ape)

source("../Just_universal/code/pairmemo.R")
pairmemo.dir = "/data-belle/Mexico_temperature/pairmemo"

library(future.apply)
plan(multiprocess)

satellite.data.dir = "data/RAW/MODIS.AQUA.TERRA.LST.NDVI/stage2"

n.folds = 10
available.years = 2003 : 2015

temp.ground.vars = c(
    "ground.temp.lo", "ground.temp.mean", "ground.temp.hi")
nontemp.ground.vars = c(
    "r.humidity.mean", "bar.mean", "wind.speed.mean", "rain.mean")

# Import the clipped grid of the Mexico study area.
fullgrid = fread("data/work/mexico_grid_ndvi_water_final.csv")
fullgrid[, lstid := paste0(long_lst, "-", lat_lst)]
fullgrid[, ndviid := paste0(long_ndvi, "-", lat_ndvi)]
fullgrid = fullgrid[,.(
    lstid, long_lst, lat_lst,
    elevation, aspectmean, roaddenmean, openplace,
    ndviid, long_ndvi, lat_ndvi,
    in_water)]
setkey(fullgrid, lstid)

# Load the data from ground stations.
ground = readRDS("data/work/all_stations_final.rds")
# Remove the `lstid` column, because it isn't formatted consistently
# with the satellite data.
ground[, lstid := NULL]
setnames(ground,
    c("low.temp", "temp.mean", "hi.temp"),
    temp.ground.vars)
# There are 39 cases in which we have more than one observation
# for a particular station and day (which all happen to be in 2013
# from station 76677). In each case, keep only the second of the
# two observations.
stopifnot(ground[, .N, by = .(date, stn)][N > 1, .N] == 39)
ground = ground[, head(.SD, 1), by = .(date, stn)]
setkey(ground, stn)
# Each station should have only one position.
stopifnot(all(
    ground[, nrow(unique(.SD)), by = stn,
       .SDcols = c("latitude", "longitude")]$V1
    == 1))

# Find the nearest `lstid` and `ndviid` (i.e., LST ID and NDVI
# ID) for each station.
nearest.id = function(longvar, latvar, idvar)
   {ground.station.pos = as.matrix(unique(ground[, .(longitude, latitude)]))
    fullgrid.pos = as.matrix(fullgrid[, c(longvar, latvar), with = F])
    neighbors = as.data.table(get.knnx(
        fullgrid.pos, ground.station.pos, k = 1))
    neighbors[, stn := unique(ground$stn)]
    neighbors[[idvar]] = fullgrid[neighbors$nn.index, idvar, with = F][[1]]
    neighbors[, c("stn", idvar), with = F]}
ground = merge(ground, by = "stn",
    nearest.id("long_lst", "lat_lst", "lstid"))
ground = merge(ground, by = "stn",
    nearest.id("long_ndvi", "lat_ndvi", "ndviid"))

get.satellite.data = function(satellite, product, the.year)
   {stopifnot(satellite %in% c("terra", "aqua"))
    stopifnot(product %in% c("temperature", "vegetation"))
    message("Loading satellite data: ", paste(satellite, product, the.year))
    d = as.data.table(readRDS(file.path(satellite.data.dir,
        sprintf("%s%s_%d.rds",
            c(terra = "MOD", aqua = "MYD")[satellite],
            c(temperature = "11A1", vegetation = "13A3")[product],
            the.year))))
    if (product == "vegetation")
        setnames(d, "lstid", "ndviid")
    message("Subsetting satellite data")
    d = (if (product == "temperature") d[
            (!is.na(d.tempc) | !is.na(n.tempc)) & lstid %in% fullgrid$lstid,
            .(lstid, yday = yday(day), d.tempc, n.tempc)]
        else d[
            !is.na(ndvi) & ndviid %in% fullgrid$ndviid,
            .(ndviid, month = month(day), ndvi)])
    message("Writing satellite data")
    d}
get.satellite.data = pairmemo(get.satellite.data, pairmemo.dir, fst = T)

model.dataset = function(the.year)
   {message("Merging data sources for ", the.year)

    message("Constructing grid")
    d = CJ(sorted = F,
        lstid = fullgrid$lstid,
        date = seq(
            as.Date(paste0(the.year, "-01-01")),
            as.Date(paste0(the.year, "-12-31")),
            by = 1))

    message("Rejiggering columns")
    d[, `:=`(
        yday = yday(date),
        month = month(date),
        date = NULL,
        ndviid = fullgrid[lstid, ndviid])]

    message("Merging in ground measurements")
    # This step can increase the number of rows in `d` a bit because
    # there can more than one ground station per `lstid`.
    d = merge(d,
        ground[
            year(date) == the.year,
            c(
                mget(c("stn", "lstid")),
                .(yday = yday(date)),
                mget(temp.ground.vars), mget(nontemp.ground.vars))],
        by = c("lstid", "yday"),
        all.x = T)

    # Merge in satellite data.
    d = with(list(),
       {for (satellite in c("terra", "aqua"))
           {st = get.satellite.data(satellite, "temperature", the.year)
            message("Merging in ", satellite, " temperature")
            d = merge(d, st, by = c("lstid", "yday"), all.x = T)
            setnames(d,
                c("d.tempc", "n.tempc"),
                paste0(satellite, c(".temp.day", ".temp.night")))
            sv = get.satellite.data(satellite, "vegetation", the.year)
            message("Merging in ", satellite, " NDVI")
            d = merge(d, sv, by = c("ndviid", "month"), all.x = T)
            setnames(d, "ndvi", paste0(satellite, ".ndvi"))}
        d})

    message("Merging in land-use data")
    d = merge(d,
        fullgrid[, .(
            lstid, elevation, aspectmean, roaddenmean, openplace)],
        by = "lstid",
        all.x = T)

    # Remove some high correlations by mean-centering.
    message("Decorrelating")
    d[, bar.mean := bar.mean - mean(bar.mean, na.rm = T), by = elevation]
    d[, roaddenmean := roaddenmean - mean(roaddenmean), by = openplace]

    # Within each grid cell, linearly interpolate missing
    # satellite temperatures on the basis of day.
    for (vname in grep("^(aqua|terra).temp\\.", colnames(d), value = T))
       {message("Interpolating ", vname)
        d[, paste0(vname, ".imputed") := is.na(get(vname))]
        d[, (vname) := approx(
                x = yday, y = get(vname), xout = yday,
                method = "linear", rule = 2)$y,
            by = lstid]}

    message("Reselecting variables")
    d = d[, .(
        lstid, yday,
        stn,
        ground.temp.lo, ground.temp.mean, ground.temp.hi,
        terra.temp.day, terra.temp.day.imputed,
        terra.temp.night, terra.temp.night.imputed,
        aqua.temp.day, aqua.temp.day.imputed,
        aqua.temp.night, aqua.temp.night.imputed,
        ndvi = (terra.ndvi + aqua.ndvi)/2,
        elevation, aspectmean, roaddenmean,
        r.humidity.mean, bar.mean, rain.mean, wind.speed.mean,
        openplace,
        time.sin = sinpi(2 * (yday - 1)/(max(yday, na.rm = T) - 1)),
        time.cos = cospi(2 * (yday - 1)/(max(yday, na.rm = T) - 1)))]

    message("Standardizing variables")
    for (col in setdiff(colnames(d), c(
            "stn", "lstid",
            temp.ground.vars,
            grep("imputed", colnames(d), value = T),
            "yday")))
        d[[col]] = arm::rescale(d[[col]])

    # Split the ground stations into cross-validation folds. This
    # has to be done on a per-year basis because not all stations
    # are present in all cases.
    #
    # The actual number of folds will be less than `n.folds` if
    # there aren't at least `n.folds` stations.
    message("Choosing cross-validation folds")
    set.seed(the.year)
    stns = sort(na.omit(unique(d$stn)))
    stn.folds = data.table(key = "stn",
        stn = stns,
        fold = sample(rep(1 : n.folds, len = length(stns))))
    d = merge(d, stn.folds, by = "stn", all = T)

    message("Setting key")
    setkey(d, stn, yday)

    message("Writing")
    d}
model.dataset = pairmemo(model.dataset, pairmemo.dir, mem = T, fst = T)

# Create a matrix `stns.by.dist` such that `stns.by.dist[lstid,]`
# gives a vector of all stations ordered by distance from
# `lstid`, with the closest first.
stns.by.dist = with(list(),
   {fullgrid.pos = as.matrix(
        fullgrid[, .(long_lst, lat_lst)])
    ground.station.pos = as.matrix(unique(
        ground[order(stn), .(longitude, latitude)]))
    neighbors = get.knnx(
        ground.station.pos, fullgrid.pos,
        k = length(unique(ground$stn)))$nn.index
    neighbors = apply(neighbors, 2, function(v) sort(unique(ground$stn))[v])
    rownames(neighbors) = fullgrid$lstid
    neighbors})

impute.nontemp.ground.vars = function(d, fold.i)
  # For each ground-station variable (other than the DV,
  # temperature), substitute values that are missing or are
  # currently in the test fold. We use the nearest station
  # that has an eligible value on the same day (or a
  # previous day, if we can't find one on the same day).
   {d = copy(d)
    for (ri in 1 : nrow(d))
       {this = d[ri,]
        for (vname in nontemp.ground.vars)
            if ((!is.null(fold.i) && this$fold == fold.i) ||
                    is.na(this[[vname]]))
               {found = F
                for (the.yday in this$yday : 1)
                   {for (the.stn in stns.by.dist[this$lstid,])
                       {other = d[.(the.stn, the.yday),]
                        if ((is.null(fold.i) || other$fold != fold.i)
                                && !is.na(other[[vname]]))
                           {d[ri, vname] = other[[vname]]
                            found = T
                            break}}
                    if (found)
                        break}
                if (!found)
                    stop("No match found for ", fold.i, " ", ri, " ", vname)}}
    d}

train.model = function(dataset)
    lmer(data = dataset, ground.temp ~
        terra.temp.day * terra.temp.day.imputed +
        terra.temp.night * terra.temp.night.imputed +
        aqua.temp.day * aqua.temp.day.imputed +
        aqua.temp.night * aqua.temp.night.imputed +
        ndvi +
        time.sin + time.cos +
        elevation + aspectmean + roaddenmean +
        r.humidity.mean + bar.mean + rain.mean + wind.speed.mean +
        openplace +
        (1 | yday))

run.cv = function(the.year, dvname)
  # Under cross-validation, predict ground temperature using
  # satellite temperature over the given year.
   {d.master = model.dataset(the.year)[!is.na(get(dvname))]
    setnames(d.master, dvname, "ground.temp")
    d.master[, setdiff(temp.ground.vars, dvname) := NULL]
    message("run.cv: ", the.year, " ", dvname)

    results = future_lapply(1 : n.folds, function(fold.i)
       {d = impute.nontemp.ground.vars(d.master, fold.i)
        m = train.model(d[fold != fold.i])
        d[fold == fold.i, .(stn, yday, pred =
            predict(m, .SD, allow.new.levels = T))]})

    for (d in results)
        d.master[.(d$stn, d$yday), pred := d$pred]

    cbind(d.master, year = the.year, dv = dvname)}
run.cv = pairmemo(run.cv, pairmemo.dir, mem = T, fst = T)

multi.run.cv = function(years)
  # Run cross-validation for each outcome in each of the given years,
  # and combine all the results into one big data.table.
   rbindlist(unlist(recursive = F, lapply(years, function(the.year)
       lapply(temp.ground.vars, function(dv)
           run.cv(the.year, dv)))))

months2seasons = factor(c(
  # From: Just, A. C., Wright, R. O., Schwartz, J., Coull, B. A.,
  # Baccarelli, A. A., Tellez-Rojo, M. M., … Kloog, I. (2015).
  # Using high-resolution satellite aerosol optical depth to
  # estimate daily PM_{2.5} geographical distribution in Mexico
  # City. Environmental Science & Technology, 49(14), 8576–8584.
  # doi:10.1021/acs.est.5b00859
    "ColdDry",  # Jan
    "ColdDry",  # Feb
    "WarmDry",  # Mar
    "WarmDry",  # Apr
    "Rainy",    # May
    "Rainy",    # Jun
    "Rainy",    # Jul
    "Rainy",    # Aug
    "Rainy",    # Sep
    "Rainy",    # Oct
    "ColdDry",  # Nov
    "ColdDry")) # Dec

summarize.cv.results = function(multirun.output)
   {d = copy(multirun.output)
    d[, dv := substr(dv, nchar("ground.temp.") + 1, 1e6)]
    d[, season := months2seasons[
        month(as.Date(paste0(year, "-01-01")) + yday)]]
    idist = 1 / as.matrix(dist(unique(
        ground[order(stn), .(latitude, longitude)])))
    diag(idist) = 0
    ustns = unique(ground$stn)
    j1 = quote(.(.N, sd = sd(ground.temp), rmse = sqrt(mean((ground.temp - pred)^2)),
        R2 = cor(ground.temp, pred)^2))
    list(
        overall = cbind(
            d
                [, eval(j1), keyby = .(year, dv)]
                [, .(year, dv,
                    N, sd, rmse, "sd - rmse" = sd - rmse, R2)],
            d
                [, .(mean.obs = mean(ground.temp), mean.pred = mean(pred)),
                    keyby = .(year, dv, stn)]
                [, .(R2.spatial = cor(mean.obs, mean.pred)^2),
                    keyby = .(year, dv)]
                [, .(R2.spatial)],
            d
                [, .(delta.obs = ground.temp - mean(ground.temp), delta.pred = pred - mean(pred)),
                    keyby = .(year, dv, stn)]
                [, .(R2.temporal = cor(delta.obs, delta.pred)^2),
                    keyby = .(year, dv)]
                [, .(R2.temporal)],
            d
                [order(year, dv, yday, stn),
                    .(p = Moran.I(
                        pred - ground.temp,
                        idist[ustns %in% stn, ustns %in% stn])$p.value),
                    by = .(year, dv, yday)]
                [, mean(p < .05), by = .(year, dv)]
                [, .("Moran ps < .05" = V1)]),
        by.season = cbind(
            d
                [, eval(j1), keyby = .(year, dv, season)]
                [, .(year, dv, season,
                    N, sd, rmse, "sd - rmse" = sd - rmse)],
            d
                [, .(stn, merr = mean(pred - ground.temp)),
                    keyby = .(year, dv, season, stn)]
                [,
                    .(Moran.I(
                        merr,
                        idist[ustns %in% stn, ustns %in% stn])$p.value),
                    by = .(year, dv, season)]
                [, .("Moran p" = V1)]))}

dedupe.lstid.days = function(d)
  # There are some `lstids` that have more than one
  # station, but we want to look up some data frames by
  # `lstid`. So just keep the first station per `lstid` and
  # day.
    d[, head(.SD, 1), by = .(lstid, yday)]

predict.temps = function(file)
  # Predict a low, mean, and high temperature for each location
  # and date.
   {orig = readRDS(file)
    d.query = with(orig, data.table(date = day, long = long48, lat = lat48))

    # Find an appropriate `lstid` for each location.
    query.pos = as.matrix(d.query[, .(long, lat)])
    fullgrid.pos = as.matrix(fullgrid[, .(long_lst, lat_lst)])
    neighbors = as.data.table(get.knnx(
        fullgrid.pos, query.pos, k = 1))
    d.query$lstid = fullgrid[neighbors$nn.index, lstid]

    f = function(slice)
       {the.year = slice[1, year(date)]

        d.model = model.dataset(the.year)
        message("Subsetting for ", the.year)
        d.model = d.model[!is.na(stn) |
            paste(lstid, yday) %in% slice[, paste(lstid, yday(date))]]
        message("Imputing non-temperature ground variables for ", the.year)
        d.model = impute.nontemp.ground.vars(d.model, fold.i = NULL)
        setkey(d.model, lstid, yday, stn)

        # If for a given `lstid` and day we have an actual ground
        # measurement, we'll use that instead of the model prediction.
        gr = ground[year(date) == the.year]
        gr[, yday := yday(date)]
        gr = dedupe.lstid.days(gr)
        setkey(gr, lstid, yday)

        sapply(simplify = F, temp.ground.vars, function(dvname)
           {message("Predicting ", dvname)
            d.model[, ground.temp := get(dvname)]
            m = train.model(d.model[!is.na(ground.temp)])
            newdata = dedupe.lstid.days(d.model)[
                .(slice$lstid, yday(slice$date))]
            pred = predict(m, newdata = newdata, allow.new.levels = T)
            observed = gr[.(slice$lstid, yday(slice$date)), get(dvname)]
            ifelse(is.na(observed), pred, observed)})}
    d.query[, (temp.ground.vars) := f(.SD),
       by = year(date),
       .SDcols = c("date", "lstid")]

    d.query[, .(
        temperature.lo = ground.temp.lo,
        temperature.mean = ground.temp.mean,
        temperature.hi = ground.temp.hi)]}
predict.temps = pairmemo(predict.temps, pairmemo.dir, mem = T, fst = T)
