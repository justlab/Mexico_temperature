# Use mixed-effects regression models to predict ground readings
# of temperature from several weather stations around Mexico City
# using satellite temperature measurements, satellite vegetation
# measurements, and some other variables. There are three DVs:
# daily high temperatures (according to the ground weather
# stations), daily lows, and daily means. A separate model is fit
# per year and DV.
#
# --- Usage ---
#
# Examine cross-validated predictive accuracy for one year and DV:
#
#     summarize.cv.results(run.cv(2012L, "ground.temp.mean"))
#
# For all DVs for a given year:
#
#     summarize.cv.results(multi.run.cv(2012L))
#
# For all DVs and all years:
#
#     summarize.cv.results(multi.run.cv(available.years))
#
# Get predictions of the lows, highs, and mean temperature for
# new positions and days in the study area:
#
#     predict.temps(FILENAME)

suppressPackageStartupMessages(
   {library(data.table)
    library(fst)
    library(FNN)
    library(lme4)
    library(optimx)
    library(ape)
    library(sf)
    library(future.apply)
    library(zeallot)
    library(caret)})

source("../Just_universal/code/pairmemo.R")
pairmemo.dir = "/data-belle/Mexico_temperature/pairmemo"

c(get.ground, in.study.area, pred.area, crs.lonlat, crs.mexico.city) %<-% local(
   {source("stations.R", local = T)
    list(get.ground, in.study.area, pred.area, crs.lonlat, crs.mexico.city)})

satellite.temperature.dir = "/data-belle/Mexico_temperature/lst_c006/mex.lst"
satellite.vegetation.dir = "/data-belle/Mexico_temperature/ndvi_c006"
elevation.path = "/data-belle/Mexico_temperature/elevation/srtm30_extracted.fst"

plan(multiprocess)

n.folds = 10
available.years = 2003 : 2017
master.grid.year = 2012L
  # This needs to be a year for which we have satellite vegetation
  # data, but the exact value shouldn't matter much.

satellite.codes = c(terra = "MOD", aqua = "MYD")

temp.ground.vars = c(
    "ground.temp.lo", "ground.temp.mean", "ground.temp.hi")
nontemp.ground.vars = c(
    "wind.speed.mean")

get.nonsatellite.data = function()
   {message("Loading master grid")
    master.grid <<- rbindlist(Map(read.vegetation.file, full.grid = T,
        grep(value = T, "\\.A\\d{4}001\\.[^/]+$",
            vegetation.paths("aqua", master.grid.year))))
    master.grid <<- master.grid[
      in.study.area(x, y),
      .(lon = x, lat = y)]
    stopifnot(nrow(unique(master.grid)) == nrow(master.grid))
    # Determine which rows of the master grid are in the prediction
    # area.
    message("Finding prediction area")
    master.grid[, in.pred.area := local(
      {mg = st_transform(crs = crs.mexico.city,
           st_as_sf(master.grid, coords = c("lon", "lat"), crs = crs.lonlat))
       pa = st_transform(crs = crs.mexico.city, pred.area())
       sti = st_intersects(mg, pa, sparse = F)
       rowSums(sti) > 0})]

    message("Loading elevation")
    # Read from fst files produced by `prepare_elevation_mex.Rmd`.
    elevation = read_fst(elevation.path, as.data.table = T)[
        in.study.area(x, y)]
    set.mrows(elevation, "x", "y")
    stopifnot(!anyDuplicated(elevation$mrow))
    master.grid[elevation$mrow, elevation := elevation$elev_filtered]
    stopifnot(!anyNA(master.grid$elevation))

    message("Loading data from ground stations")
    stations = copy(get.ground()$stations)
    ground = copy(get.ground()$obs)

    setnames(ground,
        c("temp.C.min", "temp.C.mean", "temp.C.max"),
        temp.ground.vars)
    setnames(ground, "wind.speed.mps.mean", "wind.speed.mean")
    set.mrows(stations, "lon", "lat")
    ground[, mrow := stations[.(ground$stn), mrow]]

    # Create a matrix `stns.by.dist` such that
    # `stns.by.dist[mrow,]` gives a vector of all stations
    # ordered by distance from `mrow`, with the closest first.
    message("Populating stns.by.dist")
    stns.by.dist = get.knnx(
        as.matrix(stations[order(stn), .(lon, lat)]),
        master.grid[, .(lon, lat)],
        k = nrow(stations))$nn.index
    stns.by.dist = apply(stns.by.dist, 2, function(v)
        sort(stations$stn)[v])

    list(master.grid, ground, stations, stns.by.dist)}
get.nonsatellite.data = pairmemo(get.nonsatellite.data, pairmemo.dir, mem = T)

set.mrows = function(d, longitude.col, latitude.col)
  # Adds to the given data table a column `mrow` that specifies
  # the corresponding row of `master.grid`.
   {x = d[[longitude.col]]
    y = d[[latitude.col]]
    stopifnot(all(in.study.area(x, y)))
    set(d, j = "mrow", value =
        get.knnx(master.grid[, .(lon, lat)], cbind(x, y), k = 1)$nn.index[,1])}

get.satellite.data = function(satellite, product, the.year)
   {stopifnot(satellite %in% c("terra", "aqua"))
    stopifnot(product %in% c("temperature", "vegetation"))
    message("Loading satellite data: ", paste(satellite, product, the.year))

    if (product == "temperature")
      # Read from fst files produced by
      # https://gitlab.com/ihough/modis_lst_hdf_to_fst
       {d = read_fst(
            file.path(satellite.temperature.dir,
                sprintf("%s11A1_%d.fst", satellite.codes[satellite], the.year)),
            columns = c("day", "lon", "lat", "LST_Day_1km", "LST_Night_1km"),
            as.data.table = T)
        setnames(d,
            c("LST_Day_1km", "LST_Night_1km"),
            c("temp.day", "temp.night"))
        message("Subsetting temperatures")
        d = d[(!is.na(temp.day) | !is.na(temp.night)) &
            in.study.area(lon, lat)]
        message("Setting mrows")
        set.mrows(d, "lon", "lat")
        d = d[, .(mrow, yday = yday(day), temp.day, temp.night)]}

    else
      # For vegetation, read from the original HDFs.
       {month.daynums = c(
            '001', '032', '060', '091', '121', '152',
            '182', '213', '244', '274', '305', '335')

        d = rbindlist(future_lapply(vegetation.paths(satellite, the.year), function(fpath)
           {d = read.vegetation.file(fpath)

            d = d[in.study.area(x, y)]
            set.mrows(d, "x", "y")
            d[, `:=`(x = NULL, y = NULL)]

            daynum = regmatches(fpath, regexec("\\.A\\d{4}(\\d{3})", fpath))[[1]][2]
            d$month = which(
                month.daynums == daynum |
                month.daynums == sprintf('%03d', as.integer(daynum) - 1))

            setcolorder(d, c("mrow", "month", "ndvi"))
            d}))}

    message("Writing satellite data")
    d}
get.satellite.data = pairmemo(get.satellite.data, pairmemo.dir, fst = T)

vegetation.paths = function(satellite, the.year)
    grep(fixed = T, value = T,
        paste0(satellite.codes[satellite], "13A3.A", the.year),
        dir(satellite.vegetation.dir, full.names = T))

read.vegetation.file = function(fpath, full.grid = F)
   {suppressPackageStartupMessages(
       {library(raster)
        library(rgdal)
        library(gdalUtils)})
    subdataset = paste0(
        "HDF4_EOS:EOS_GRID:",
        fpath,
        ":MOD_Grid_monthly_1km_VI:1 km monthly NDVI")
    d = raster(readGDAL(subdataset, silent = T))
    if (full.grid)
        d[,] = 1
    d = as.data.table(spTransform(rasterToPoints(d, spatial = T),
        "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
    if (!full.grid)
       {setnames(d, "band1", "ndvi")
        # The scale factor has already been applied, but by
        # multiplication instead of division, so divide by the
        # square.
        gi = paste(gdalinfo(subdataset), collapse = " ")
        scale.factor = as.numeric(regmatches(gi,
            regexec(" scale_factor=(\\d+)", gi))[[1]][2])
        d$ndvi = d$ndvi / scale.factor^2
        # Missing points have already been thrown out.
        stopifnot(!anyNA(d))}
    d}

mrow.sets = list()

model.dataset = function(the.year, mrow.set = NULL, nonmissing.ground.temp = F)
   {if (!xor(!is.null(mrow.set), nonmissing.ground.temp))
        stop("Either enable `nonmissing.ground.temp` to get all points with non-missing ground temperatures, or choose an `mrow.set` to get all days for the selected `mrow`s.")

    message("Merging data sources: ",
        the.year, " ", mrow.set, " ", nonmissing.ground.temp)

    message("Constructing grid")
    d = CJ(sorted = F,
        mrow = (if (nonmissing.ground.temp)
                sort(ground[year(date) == the.year, unique(mrow)])
            else
                mrow.sets[[mrow.set]]),
        date = seq(
            as.Date(paste0(the.year, "-01-01")),
            as.Date(paste0(the.year, "-12-31")),
            by = 1))

    message("Rejiggering columns")
    d[, `:=`(
        yday = yday(date),
        month = month(date),
        date = NULL,
        elevation = master.grid[mrow, elevation])]

    message("Merging in ground measurements")
    # This step can increase the number of rows in `d` a bit because
    # there can more than one ground station per `mrow`.
    d = merge(d,
        ground[
            year(date) == the.year,
            c(
                mget(c("stn", "mrow")),
                .(yday = yday(date)),
                mget(temp.ground.vars), mget(nontemp.ground.vars))],
        by = c("mrow", "yday"),
        all.x = T)

    # Merge in satellite data.
    d = local(
       {for (satellite in c("terra", "aqua"))
           {st = get.satellite.data(satellite, "temperature", the.year)
            message("Merging in ", satellite, " temperature")
            d = merge(d, st, by = c("mrow", "yday"), all.x = T)
            setnames(d,
                c("temp.day", "temp.night"),
                paste0(satellite, c(".temp.day", ".temp.night")))
            sv = get.satellite.data(satellite, "vegetation", the.year)
            message("Merging in ", satellite, " NDVI")
            d = merge(d, sv, by = c("mrow", "month"), all.x = T)
            # Very rarely, NDVI is missing. In such cases, set it
            # to the mean of nonmissing values on that month.
            if (anyNA(d$ndvi))
               {message("Imputing ", sum(is.na(d$ndvi)), " missing values")
                d[, by = month, ndvi := ifelse(is.na(ndvi),
                     mean(ndvi, na.rm = T),
                     ndvi)]}
            setnames(d, "ndvi", paste0(satellite, ".ndvi"))}
        d})

    combine.sat = function(terra, aqua)
        ifelse(is.na(terra),
            ifelse(is.na(aqua),
                NA,
                aqua),
            ifelse(is.na(aqua),
                terra,
                (terra + aqua)/2))
    d[, `:=`(
        satellite.temp.day = combine.sat(terra.temp.day, aqua.temp.day),
        satellite.temp.night = combine.sat(terra.temp.night, aqua.temp.night))]

    # Within each grid cell, linearly interpolate missing
    # satellite temperatures on the basis of day.
    for (vname in c("satellite.temp.day", "satellite.temp.night"))
       {message("Interpolating ", vname)
        d[, paste0(vname, ".imputed") := is.na(get(vname))]
        d[, (vname) := approx(
                x = yday, y = get(vname), xout = yday,
                method = "linear", rule = 2)$y,
            by = mrow]}

    if (nonmissing.ground.temp)
       {message("Dropping")
        d = d[0 < rowSums(!is.na(d[, ..temp.ground.vars]))]}

    message("Reselecting variables")
    max.yday = 365 + lubridate::leap_year(the.year)
    d = d[, .(
        mrow, yday,
        stn,
        ground.temp.lo, ground.temp.mean, ground.temp.hi,
        satellite.temp.day, satellite.temp.day.imputed,
        satellite.temp.night, satellite.temp.night.imputed,
        ndvi = (terra.ndvi + aqua.ndvi)/2,
        elevation,
        wind.speed.mean,
        time.sin = sinpi(2 * (yday - 1)/(max.yday - 1)),
        time.cos = cospi(2 * (yday - 1)/(max.yday - 1)))]

    if (nonmissing.ground.temp)
      # Split the ground stations into cross-validation folds. This
      # has to be done on a per-year basis because not all stations
      # are present in all cases. Only stations that are in the
      # prediction area get a fold; other stations are always
      # used for training (their fold is set to -1).
      #
      # The actual number of folds will be less than `n.folds` if
      # there aren't at least `n.folds` stations.
       {message("Choosing cross-validation folds")
        set.seed(the.year)
        stns = setdiff(sort(na.omit(unique(d$stn))),
            stations[master.grid[stations$mrow, !in.pred.area], stn])
        stn.folds = data.table(key = "stn",
            stn = stns,
            fold = sample(rep(1 : n.folds, len = length(stns))))
        d = merge(d, stn.folds, by = "stn", all = T)
        d[is.na(fold), fold := -1]}
    else
        d$fold = NA

    message("Setting key")
    setkey(d, stn, yday)

    message("Writing")
    d}
model.dataset = pairmemo(model.dataset, pairmemo.dir, mem = T, fst = T)

impute.nontemp.ground.vars = function(d.orig, fold.i)
  # For each ground-station variable (other than the DV,
  # temperature), substitute values that are missing or are
  # currently in the test fold. We use the nearest station
  # that has an eligible value on the same day (or a
  # previous day, if we can't find one on the same day).
   {d = copy(d.orig)
    folds = d$fold
    ydays = d$yday
    mrows = d$mrow
    ustn = as.character(unique(d.orig$stn))

    for (vname in nontemp.ground.vars)
       {column = d[[vname]]

        # Create `other.vs` (and likewise `other.folds`) such that
        # `other.vs[[stn]][yday]` gives the value of `vname` at
        # station `stn` on day of the year `yday`.
        other.vs = new.env(parent = emptyenv())
        other.folds = new.env(parent = emptyenv())
        for (the.stn in ustn)
           {piece = d.orig[
                stn == as.integer(the.stn),
                .(yday, v = get(vname), fold)]
            other.vs[[the.stn]] = rep(NA_real_, max(d$yday))
            other.vs[[the.stn]][piece$yday] = piece$v
            other.folds[[the.stn]] = rep(NA_integer_, max(d$yday))
            other.folds[[the.stn]][piece$yday] = piece$fold}

        # Now look at each value of `column` in turn to see if it
        # needs a substitute.
        for (ri in seq_along(column))
           {if ((!is.null(fold.i) && folds[ri] == fold.i) ||
                    is.na(column[ri]))
               {found = F
                for (the.yday in ydays[ri] : 1)
                   {for (the.stn in as.character(stns.by.dist[mrows[ri],]))
                       {if (exists(the.stn, other.vs)
                                && (is.null(fold.i) || other.folds[[the.stn]][the.yday] != fold.i)
                                && !is.na(other.vs[[the.stn]][the.yday]))
                           {column[ri] = other.vs[[the.stn]][the.yday]
                            found = T
                            break}}
                    if (found)
                        break}
                if (!found)
                    stop("No match found for ", fold.i, " ", ri, " ", vname)}}

        d[[vname]] = column}

    d}

train.model = function(dataset)
   {fe = (~
        satellite.temp.day + satellite.temp.day.imputed +
        satellite.temp.night + satellite.temp.night.imputed +
        ndvi +
        time.sin + time.cos +
        elevation +
        wind.speed.mean)
    preproc = preProcess(method = c("center", "scale"),
        dataset[,
            setdiff(all.vars(fe), grep("imputed", all.vars(fe), val = T)),
            with = F])
    m = lmer.alternatives(
        data = predict(preproc, dataset),
        formula = update.formula(fe, ground.temp ~ . +
            (1 + satellite.temp.day + satellite.temp.night | yday)),
        optimizers = list(
            list(),
            list(control = lmerControl(optimizer = "Nelder_Mead")),
            list(control = lmerControl(optimizer = "optimx",
                optCtrl = list(method = "L-BFGS-B"))),
            list(control = lmerControl(optimizer = "optimx",
                optCtrl = list(method = "nlminb")))))
    function(newdata)
       predict(m, newdata = predict(preproc, newdata),
           allow.new.levels = T)}

lmer.alternatives = function(formula, data, optimizers)
  # Try calling `lmer` with each list of arguments in `optimizers`
  # and return the model from the first fit that doesn't produce any
  # warnings.
   {for (arglist in optimizers)
        tryCatch(
            return(do.call(lmer, c(list(formula, data), arglist))),
            warning = function(w) NULL)
    stop("All lmer.alternatives failed")}

run.cv = function(the.year, dvname)
  # Under cross-validation, predict ground temperature using
  # satellite temperature over the given year.
   {d.master = copy(model.dataset(the.year, nonmissing.ground.temp = T))
    setnames(d.master, dvname, "ground.temp")
    d.master[, setdiff(temp.ground.vars, dvname) := NULL]
    message("run.cv: ", the.year, " ", dvname)

    bar = txtProgressBar(min = 0, max = n.folds, style = 3)
    for (fold.i in 1 : n.folds)
       {d = impute.nontemp.ground.vars(d.master, fold.i)
        f.pred = train.model(d[fold != fold.i])
        d.master[fold == fold.i, pred := f.pred(d[fold == fold.i])]
        setTxtProgressBar(bar, fold.i)}
    close(bar)

    cbind(d.master, year = the.year, dv = dvname)}
run.cv = pairmemo(run.cv, pairmemo.dir, mem = T, fst = T)

multi.run.cv = function(years)
  # Run cross-validation for each outcome in each of the given years,
  # and combine all the results into one big data.table.
  {args = expand.grid(the.year = years, dv = temp.ground.vars,
       stringsAsFactors = F)
   rbindlist(lapply(1 : nrow(args), function(i)
       run.cv(args[i, "the.year"], args[i, "dv"])))}

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
   {d = multirun.output[fold != -1]
    d[, dv := substr(dv, nchar("ground.temp.") + 1, 1e6)]
    d[, season := months2seasons[
        month(as.Date(paste0(year, "-01-01")) + yday)]]
    idist = 1 / as.matrix(dist(stations[, .(lon, lat)]))
    diag(idist) = 0
    ustns = stations$stn
    j1 = quote(.(.N, stn = length(unique(stn)),
        sd = sd(ground.temp), rmse = sqrt(mean((ground.temp - pred)^2)),
        R2 = cor(ground.temp, pred)^2))
    list(
        overall = cbind(
            d
                [, eval(j1), keyby = .(year, dv)]
                [, .(year, dv,
                    N, stn, sd, rmse, "sd - rmse" = sd - rmse, R2)],
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
        by.imp =
            d
                [, eval(j1), keyby = .(year, dv,
                    imp.d = satellite.temp.day.imputed,
                    imp.n = satellite.temp.night.imputed)]
                [, .(year, dv, imp.d, imp.n,
                    N, stn, sd, rmse, "sd - rmse" = sd - rmse)],
        by.season = cbind(
            d
                [, eval(j1), keyby = .(year, dv, season)]
                [, .(year, dv, season,
                    N, stn, sd, rmse, "sd - rmse" = sd - rmse)],
            d
                [, .(stn, merr = mean(pred - ground.temp)),
                    keyby = .(year, dv, season, stn)]
                [,
                    .(Moran.I(
                        merr,
                        idist[ustns %in% stn, ustns %in% stn])$p.value),
                    by = .(year, dv, season)]
                [, .("Moran p" = V1)]))}

dedupe.mrow.days = function(d)
  # There are some `mrows` that have more than one
  # station, but we want to look up some data frames by
  # `mrow`. So just keep the first station per `mrow` and
  # day.
    d[, head(.SD, 1), by = .(mrow, yday)]

predict.temps = function(file)
  # Predict a low, mean, and high temperature for each location
  # and date.
   {orig = readRDS(file)
    d.query = with(orig, data.table(date = day, lon = long48, lat = lat48))

    message("Setting mrows")
    set.mrows(d.query, "lon", "lat")

    f = function(slice)
       {the.year = slice[1, year(date)]

        # Construct `d.model` to contain to have a row for each time
        # (`yday`) and place (`mrow`) that either:
        # - has a non-missing ground temperature (used for training)
        # - is present in `slice` (used for prediction)
        mrow.sets[[paste(file, the.year)]] <<- sort(unique(slice$mrow))
        d.model = model.dataset(the.year, nonmissing.ground.temp = T)
        d.model = rbind(d.model,
            model.dataset(the.year, mrow.set = paste(file, the.year))[
                paste(mrow, yday) %in% slice[, paste(mrow, yday(date))] &
                !(paste(mrow, yday) %in% d.model[, paste(mrow, yday)])])
        d.model[, ground.temp := NA_real_]
        setkey(d.model, stn, yday)
        message("Imputing non-temperature ground variables for ", the.year)
        d.model = impute.nontemp.ground.vars(d.model, fold.i = NULL)
        setkey(d.model, mrow, yday, stn)

        # If for a given `mrow` and day we have an actual ground
        # measurement, we'll use that instead of the model prediction.
        gr = ground[year(date) == the.year]
        gr[, yday := yday(date)]
        gr = dedupe.mrow.days(gr)
        setkey(gr, mrow, yday)

        sapply(simplify = F, temp.ground.vars, function(dvname)
           {message("Predicting ", dvname)
            d.model[, ground.temp := get(dvname)]
            f.pred = train.model(d.model[!is.na(ground.temp)])
            pred = f.pred(dedupe.mrow.days(d.model)[
                .(slice$mrow, yday(slice$date))])
            observed = gr[.(slice$mrow, yday(slice$date)), get(dvname)]
            ifelse(is.na(observed), pred, observed)})}
    d.query[, (temp.ground.vars) := f(.SD),
       by = year(date),
       .SDcols = c("date", "mrow")]

    d.query[, .(
        temperature.lo = ground.temp.lo,
        temperature.mean = ground.temp.mean,
        temperature.hi = ground.temp.hi)]}
predict.temps = pairmemo(predict.temps, pairmemo.dir, mem = T, fst = T)

c(master.grid, ground, stations, stns.by.dist) %<-% get.nonsatellite.data()
