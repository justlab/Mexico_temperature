library(data.table)
library(FNN)
library(lme4)
library(ape)

source("../Just_universal/code/pairmemo.R")
pairmemo.dir = "/data-belle/Mexico_temperature/pairmemo"

library(future.apply)
plan(multiprocess)

aqua.dir = "data/RAW/MODIS.AQUA.TERRA.LST.NDVI/stage2"

n.folds = 10
available.years = 2003 : 2015

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

# Load the data from ground stations.
ground = readRDS("data/work/all_stations_final.rds")
# There are 39 cases in which we have more than one observation
# for a particular station and day (which all happen to be in 2013
# from station 76677). In each case, keep only the second of the
# two observations.
stopifnot(ground[, .N, by = .(date, stn)][N > 1, .N] == 39)
ground = ground[, head(.SD, 1), by = .(date, stn)]
setkey(ground, stn)

# Find the nearest LST and NDVI ID for each station.
nearest.id = function(longvar, latvar, idvar)
   {ground.station.pos = as.matrix(unique(ground[, .(longitude, latitude)]))
    fullgrid.pos = as.matrix(fullgrid[, c(longvar, latvar), with = F])
    neighbors = as.data.table(get.knnx(
        fullgrid.pos, ground.station.pos, k = 1))
    neighbors[, stn := unique(ground$stn)]
    neighbors[[idvar]] = fullgrid[neighbors$nn.index, idvar, with = F][[1]]
    neighbors = neighbors[, c("stn", idvar), with = F]
    setkey(neighbors, stn)
    neighbors}
nearest.ids = merge(
    nearest.id("long_lst", "lat_lst", "lstid"),
    nearest.id("long_ndvi", "lat_ndvi", "ndviid"))

aqua.temp = function(year)
   {message("Loading Aqua temperature ", year)
    aqua = readRDS(file.path(aqua.dir, sprintf("MYD11A1_%d.rds", year)))
    message("Subsetting")
    d = as.data.table(aqua[
        aqua$lstid %in% nearest.ids$lstid,
        c("lstid", "day", "d.tempc", "n.tempc")])
    d = d[!is.na(d.tempc) | !is.na(n.tempc),
       .(lstid, yday = yday(day), d.tempc, n.tempc)]
    message("Writing")
    d}
aqua.temp = pairmemo(aqua.temp, pairmemo.dir, mem = T, fst = T)

aqua.ndvi = function(the.year)
   {message("Loading Aqua NDVI ", the.year)
    aqua = readRDS(file.path(aqua.dir, sprintf("MYD13A3_%d.rds", the.year)))
    colnames(aqua)[colnames(aqua) == "lstid"] = "ndviid"
    message("Subsetting")
    d = as.data.table(aqua[
        aqua$ndviid %in% nearest.ids$ndviid,
        c("ndviid", "day", "ndvi")])
    d = d[!is.na(ndvi), .(ndviid, month = month(day), ndvi)]
    message("Writing")
    d}
aqua.ndvi = pairmemo(aqua.ndvi, pairmemo.dir, mem = T, fst = T)

model.dataset = function(the.year, ground.temp.var, satellite.temp.var)
   {# Get ground temperature.
    d = ground[
        year(date) == the.year,
        c(ground.temp.var, "stn", "date", nontemp.ground.vars),
        with = F]
    setnames(d, ground.temp.var, "ground.temp")
    d[, yday := yday(date)]
    d[, month := month(date)]

    # Merge in land-use variables.
    d = merge(d, nearest.ids, by = "stn", all.x = T)
    d = merge(d,
        fullgrid[, .(
            lstid, elevation, aspectmean, roaddenmean, openplace)],
        by = "lstid",
        all.x = T)

    # Merge in satellite data.
    d = merge(d,
        aqua.temp(the.year)[
            lstid %in% d$lstid,
            c("lstid", "yday", satellite.temp.var),
            with = F],
        by = c("lstid", "yday"),
        all = T)
    setnames(d, satellite.temp.var, "satellite.temp")
    d = merge(d,
        aqua.ndvi(the.year),
        by = c("ndviid", "month"),
        all.x = T)

    # Within each grid cell, linearly interpolate missing
    # satellite temperatures on the basis of day.
    d[, satellite.temp.imputed := is.na(satellite.temp)]
    d[,
        satellite.temp := approx(
            x = yday, y = satellite.temp, xout = yday,
            method = "linear", rule = 2)$y,
        by = lstid]

    # We can now throw out rows with missing ground temperatures,
    # and keep only the columns we want.
    d = d[!is.na(ground.temp), .(
        stn, satellite.temp.imputed,
        ground.temp, satellite.temp, ndvi,
        elevation, aspectmean, roaddenmean,
        r.humidity.mean, bar.mean, rain.mean, wind.speed.mean,
        openplace, yday)]

    # Standardize most variables.
    for (col in setdiff(colnames(d), c("lstid",
            "ground.temp", "stn", "satellite.temp.imputed", "yday")))
        d[[col]] = arm::rescale(d[[col]])

    # Split the ground stations into cross-validation folds. This
    # has to be done on a per-year, per-satellite-temperature
    # basis because not all stations are present in all cases.
    #
    # The actual number of folds will be less than `n.folds` if
    # there aren't at least `n.folds` stations.
    set.seed(the.year * 10 + (satellite.temp.var == "n.tempc"))
    stns = sort(unique(d$stn))
    stn.folds = data.table(key = "stn",
        stn = stns,
        fold = sample(rep(1 : n.folds, len = length(stns))))
    d = merge(d, stn.folds, by = "stn")
    setkey(d, stn, yday)

    d}

train.model = function(dataset)
    lmer(data = dataset, ground.temp ~
        satellite.temp + ndvi +
        elevation + aspectmean + roaddenmean +
        r.humidity.mean + bar.mean + rain.mean + wind.speed.mean +
        openplace +
        (1 + satellite.temp | yday))

# Create a named list `other.stns.by.dist` such that
# `other.stns.by.dist[[stn]]` gives a vector of all the other
# stations ordered by distance from ``stn``, with the closest
# first.
stn.dists = as.matrix(dist(ground[
    , head(.SD, 1), by = stn][
    order(stn), .(latitude, longitude)]))
other.stns.by.dist = lapply(unique(ground$stn), function(the.stn)
    unique(ground$stn)[order(
       stn.dists[unique(ground$stn) == the.stn,])][-1])
names(other.stns.by.dist) = unique(ground$stn)

run.cv = function(the.year, ground.temp.var, satellite.temp.var)
  # Under cross-validation, predict ground temperature using
  # satellite temperature over the given year.
   {d.master = model.dataset(the.year, ground.temp.var, satellite.temp.var)
    message(the.year, " ", ground.temp.var, " ", satellite.temp.var)

    results = future_lapply(1 : n.folds, function(fold.i)
       {d = copy(d.master)

        # For each ground-station variable (other than the DV,
        # temperature), substitute values that are missing or are
        # currently in the test fold. We use the nearest station
        # that has an eligible value on the same day (or a
        # previous day, if we can't find one on the same day).
        for (ri in 1 : nrow(d))
           {this = d[ri,]
            for (vname in nontemp.ground.vars)
                if (this$fold == fold.i || is.na(this[[vname]]))
                   {found = F
                    for (the.yday in this$yday : 1)
                       {for (other.stn in other.stns.by.dist[[as.character(this$stn)]])
                           {other = d[.(other.stn, the.yday),]
                            if (other$fold != fold.i && !is.na(other[[vname]]))
                               {d[ri, vname] = other[[vname]]
                                found = T
                                break}}
                        if (found)
                            break}
                    stopifnot(found)}}

        m = train.model(d[fold != fold.i])
        d[fold == fold.i, .(stn, yday, pred =
            predict(m, .SD, allow.new.levels = T))]})

    for (d in results)
        d.master[.(d$stn, d$yday), pred := d$pred]

    d.master}
run.cv = pairmemo(run.cv, pairmemo.dir, mem = T, fst = T)

multirun = function(years)
  # Run cross-validation for each outcome in each of the given years,
  # and combine all the results into one big data.table.
   rbindlist(unlist(recursive = F, lapply(years, function(the.year) list(
       cbind(run.cv(the.year, "low.temp", "n.tempc"), year = the.year, dv = "lo"),
       cbind(run.cv(the.year, "temp.mean", "n.tempc"), year = the.year, dv = "mean"),
       cbind(run.cv(the.year, "hi.temp", "d.tempc"), year = the.year, dv = "hi")))))

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

summarize.results = function(multirun.output)
   {d = copy(multirun.output)
    d[, season := months2seasons[
        month(as.Date(paste0(year, "-01-01")) + yday)]]
    idist = 1 / stn.dists
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
        by.imp = d
            [, eval(j1), keyby = .(year, dv, satellite.temp.imputed)]
            [, .(year, dv, imp = satellite.temp.imputed,
                N, sd, rmse, "sd - rmse" = sd - rmse)],
        by.season = d
            [, eval(j1), keyby = .(year, dv, season)]
            [, .(year, dv, season,
                N, sd, rmse, "sd - rmse" = sd - rmse)])}
