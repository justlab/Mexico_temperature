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
# all days and squares in the prediction area:
#
#     predict.temps(2012L, "pred.area")

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
    library(caret)
    library(stringr)})

source("common.R")

satellite.temperature.dir = "/data-belle/Mexico_temperature/lst_c006/mex.lst"
satellite.vegetation.dir = "/data-belle/Mexico_temperature/ndvi_c006"
elevation.path = "/data-belle/Mexico_temperature/elevation/srtm30_extracted.fst"
mexico.city.agebs.path = "~/Jdrive/PM/Just_Lab/projects/airmex/data/gis/gisdata/AGEBS_CDMX_2010.shp"
all.agebs.path = "/data-belle/Mexico_temperature/agebs_2010"
all.agebs.year = 2010L
population.path.fmt = "~/Jdrive/PM/Just_Lab/projects/airmex/data/population/%s_cuadratic_csv"

plan(multiprocess)

n.folds = 10
available.years = 2003 : 2018
master.grid.year = 2012L
  # This needs to be a year for which we have satellite vegetation
  # data, but the exact value shouldn't matter much.

satellite.codes = c(terra = "MOD", aqua = "MYD")
satellite.tiles = c("h08v06", "h08v07")
  # https://modis-land.gsfc.nasa.gov/MODLAND_grid.html

temp.ground.vars = c(
    "ground.temp.lo", "ground.temp.mean", "ground.temp.hi")
nontemp.ground.vars = c(
    "wind.speed.mean")

get.nonsatellite.data = function()
   {message("Loading master grid")
    master.grid <<- rbindlist(Map(read.vegetation.file, full.grid = T,
        grep(value = T, "\\.A\\d{4}001\\.[^/]+$",
            vegetation.paths("aqua", master.grid.year))))
    master.grid <<- master.grid[in.study.area(lon, lat)]
    stopifnot(nrow(unique(master.grid[, .(lon, lat)])) == nrow(master.grid))
    master.grid[, mrow := .I]
    setcolorder(master.grid, "mrow")
    # Determine which region each cell of the master grid is in.
    message("Finding regions")
    master.grid[, region := local(
      {mg = st_transform(crs = crs.mexico.city,
           st_as_sf(master.grid, coords = c("lon", "lat"), crs = crs.lonlat))
       pa = st_transform(crs = crs.mexico.city, pred.area())
       sti = st_intersects(mg, pa, sparse = F)
       # Each cell of the master grid should be in at most one region.
       stopifnot(all(rowSums(sti) %in% c(0, 1)))
       factor(labels = pa$NOMGEO, apply(sti, 1, function(v)
            if (any(v)) which(v) else NA))})]
    master.grid[, in.pred.area := !is.na(region)]

    message("Loading elevation")
    # Read from fst files produced by `prepare_elevation_mex.Rmd`.
    elevation = read_fst(elevation.path, as.data.table = T)[
        in.study.area(x, y)]
    set.mrows(elevation, "x", "y")
    stopifnot(!anyDuplicated(elevation$mrow))
    master.grid[elevation$mrow, elevation := elevation$elev_filtered]
    stopifnot(!anyNA(master.grid$elevation))

    message("Loading data from ground stations")
    stations = as.data.table(get.ground()$stations)
    ground = as.data.table(get.ground()$obs)

    setkey(stations, stn)
    stations[, network := factor(network)]

    ground[, date := as.Date(date)]
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
      # The dates are actually in UTC, not our desired working
      # time zone of UTC-06:00. However, I checked a year's worth
      # of overpass times and it seems to work out that each
      # overpass is assigned to the correct UTC-06:00 date.
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
   {paths = function()
        str_subset(
            list.files(satellite.vegetation.dir, full.names = T),
            fixed(paste0(satellite.codes[satellite], "13A3.A", the.year)))
    # Check that there's at least one file for this year. If there
    # isn't, download all the files we need for this year.
    if (!length(paths()))
        download.vegetation(satellite, the.year)
    paths()}

download.vegetation = function(satellite, the.year)
   {message("Downloading vegetation for ", satellite, " ", the.year)
    suppressPackageStartupMessages(library(httr))

    creds = Sys.getenv(names = F,
        c("EARTHDATA_USERNAME", "EARTHDATA_PASSWORD"))
    if (any(creds == ""))
        stop("You need to set the environment variables EARTHDATA_USERNAME and EARTHDATA_PASSWORD. If you don't have an account, you can get one at https://urs.earthdata.nasa.gov/users/new")

    base.url = sprintf("https://e4ftl01.cr.usgs.gov/MOL%s/%s13A3.006",
        substr(toupper(satellite), 1, 1),
        satellite.codes[satellite])

    for (the.month in 1 : 12)
       {the.dir = sprintf("%s/%d.%02d.01",
            base.url, the.year, the.month)
        page = GET(the.dir)
        stop_for_status(page)

        for (tile in satellite.tiles)
           {fname = str_match(content(page, "text"),
                sprintf('<a href="([^"]+?\\.%s\\.[^"]+\\.hdf)"', tile))[,2]
            message("Getting ", fname)
            r = GET(paste0(the.dir, "/", fname),
                authenticate(creds[1], creds[2]))
            stop_for_status(r)
            writeBin(content(r, "raw"),
                file.path(satellite.vegetation.dir, fname))}}}

read.vegetation.file = function(fpath, full.grid = F)
   {suppressPackageStartupMessages(
       {library(raster)
        library(rgdal)
        library(gdalUtils)})
    subdataset = paste0(
        "HDF4_EOS:EOS_GRID:",
        fpath,
        ":MOD_Grid_monthly_1km_VI:1 km monthly NDVI")
    g = readGDAL(subdataset, silent = T)
    if (full.grid)
        g$band1 = 1
    g = as(g, "SpatialPointsDataFrame")
    d = as.data.table(spTransform(g,
        "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
    if (full.grid)
       {d[, band1 := NULL]
        setnames(d, c("x", "y"), c("lon", "lat"))
        d[, c("x_sinu", "y_sinu") := as.data.frame(g)[, c("x", "y")]]
        d = d[order(lon, lat)]}
    else
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
    withCallingHandlers(
       {for (vname in c("satellite.temp.day", "satellite.temp.night"))
           {message("Interpolating ", vname)
            d[, paste0(vname, ".imputed") := is.na(get(vname))]
            d[, (vname) := approx(
                    x = yday, y = get(vname), xout = yday,
                    method = "linear", rule = 2)$y,
                by = mrow]}},
        # Ignore warnings that `yday` is not unique within `mrow`.
        warning = function(w)
            if (grepl("collapsing to unique", conditionMessage(w)))
                invokeRestart("muffleWarning"))

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
        season = factor(months2seasons[month]),
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

impute.nontemp.ground.vars = function(d.orig, fold.i, progress = F)
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
       {column = d[, get(vname)]

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
        if (progress)
            bar = txtProgressBar(min = 0, max = length(column), style = 3)
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
                    stop("No match found for ", fold.i, " ", ri, " ", vname)}
            if (progress) setTxtProgressBar(bar, ri)}
        if (progress) close(bar)

        d[, (vname) := column]}

    d}

train.model = function(dataset)
   {fe = (~
        season * (satellite.temp.day + satellite.temp.night) +
        satellite.temp.day.imputed +
        satellite.temp.night.imputed +
        ndvi +
        time.sin + time.cos +
        elevation +
        wind.speed.mean)
    preproc = preProcess(method = c("center", "scale"),
        dataset[, mget(grep(invert = T, val = T,
            "imputed|season", all.vars(fe)))])
    m = lmer(
        data = predict(preproc, dataset),
        formula = update.formula(fe, ground.temp ~ . +
            (1 | yday)),
        control = lmerControl(check.conv.singular = "warning"))
    function(newdata)
       predict(m, newdata = predict(preproc, newdata),
           allow.new.levels = T)}

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
    idist = 1 / as.matrix(dist(stations[, .(lon, lat)]))
    diag(idist) = 0
    ustns = stations$stn
    f.r2 = function(pred, obs)
        1 - mean((obs - pred)^2)/var(obs)
    j1 = quote(.(.N, stn = length(unique(stn)),
        sd = sd(ground.temp), rmse = sqrt(mean((ground.temp - pred)^2)),
        R2 = f.r2(pred, ground.temp)))
    d[, c("lon", "lat") := stations[d$stn, .(lon, lat)]]
    d[, region := master.grid[stations[d$stn, mrow], region]]
    d[, network := stations[d$stn, network]]
    list(
        overall = cbind(
            d
                [, eval(j1), keyby = .(year, dv)]
                [, .(year, dv,
                    N, stn, sd, rmse, R2)],
            d
                [, .(ground.temp, pred, lon, lat, yday,
                        m = mean(ground.temp)),
                    by = .(year, dv)]
                [, .(var.s = mean((ground.temp - m)^2),
                        mse.s = mean((ground.temp - pred)^2)),
                    by = .(year, dv, yday, x = floor(4 * lon), y = floor(4 * lat))]
                [, .(sd.s = sqrt(mean(var.s)), rmse.s = sqrt(mean(mse.s))),
                    keyby = .(year, dv)]
                [, .(sd.s, rmse.s)],
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
                [, .("Moran p" = V1)]),
        by.region =
            d
                [year == 2018, eval(j1), keyby = .(dv, region)]
                [, .(dv, region,
                    N, stn, sd, rmse, "sd - rmse" = sd - rmse)],
        by.network =
            d
                [year == 2018, eval(j1), keyby = .(dv, network)]
                [, .(dv, network,
                    N, stn, sd, rmse, "sd - rmse" = sd - rmse)])}

predict.temps = function(the.year, mrow.set)
  # Predict a low, mean, and high temperature for every cell in
  # `mrow.set` and day in `the.year`.
   {d = model.dataset(the.year, nonmissing.ground.temp = T)
    d = rbind(
        d,
        model.dataset(the.year, mrow.set = mrow.set)[
            !(paste(mrow, yday) %in% d[, paste(mrow, yday)])])
    d[, ground.temp := NA_real_]

    message("Imputing non-temperature ground variables")
    d = impute.nontemp.ground.vars(d, fold.i = NULL, progress = T)

    for (dvname in temp.ground.vars)
       {message("Training for ", dvname)
        d[, ground.temp := get(dvname)]
        f.pred = train.model(d[!is.na(ground.temp)])
        message("Predicting ", dvname)
        d[, paste0("pred.", dvname) := f.pred(.SD)]}

    d = d[mrow %in% mrow.sets[[mrow.set]],
        keyby = .(mrow, yday),
        .SDcols = paste0("pred.", temp.ground.vars),
        head(.SD, 1)]}
predict.temps = pairmemo(predict.temps, pairmemo.dir, fst = T)

predict.temps.at = function(fname, date.col, lon.col, lat.col)
   {d = as.data.table(readRDS(fname)[,
        c(date.col, lon.col, lat.col)])
    setnames(d, c("date", "lon", "lat"))
    set.mrows(d, "lon", "lat")
    d[, yday := yday(date)]
    message("Getting model predictions")
    d[, by = .(y = year(date)), paste0("model.", temp.ground.vars) :=
        predict.temps(y, "pred.area")[.(.SD$mrow, .SD$yday),
            mget(paste0("pred.", temp.ground.vars))]]

    # Also find the temperature statistics of the nearest ground
    # station.
    message("Getting nearest temperatures")
    g = copy(ground)
    setkey(g, stn, date)
    d[, paste0("nearest.", temp.ground.vars) :=
        rbindlist(lapply(1 : .N, function(r)
           {the.date = date[r]
            for (the.stn in stns.by.dist[mrow[r],])
               {gr = g[.(the.stn, the.date), mget(temp.ground.vars)]
                if (!anyNA(gr))
                    return(gr)}
            stop()}))]

    # Also get some crude predictions of precipitation, via
    # inverse-distance weighting (IDW).
    message("Predicting precipitation")
    precip = merge(
        get.ground()$obs[, .(stn, date, precipitation.mm.total)],
        get.ground()$stations[, .(stn, lon, lat)],
        all.x = T, all.y = F,
        by = "stn")[!is.na(precipitation.mm.total)]
    d[, by = date, precipitation.mm := gstat::idw(
        formula = precipitation.mm.total ~ 1,
        locations = ~ lon + lat,
        data = precip[date == .BY$date],
        newdata = .SD,
        debug.level = 0)[, "var1.pred"]]

    d[, mget(c(
        paste0("model.", temp.ground.vars),
        paste0("nearest.", temp.ground.vars),
        "precipitation.mm"))]}

deleg.weighted.preds = function()
  # Predict temperature per day and delegación (subregion of
  # Mexico City), weighted by the population of each delegación.
  # Predictions are made with weights from people of all ages and
  # from people aged 65 and above.
   {message("Loading population counts")
    d = rbindlist(unlist(rec = F, lapply(c(F, T), function(only.65plus) lapply(
        list.files(full.names = T, sprintf(population.path.fmt,
            (if (only.65plus) "65" else "Total"))),
        function(s)
           {d = subset(fread(s), select = -TOTAL)
            d[, date := as.Date(date, format = "%m/%d/%Y")]
            d = d[year(date) %in% available.years]
            setnames(d, names(d), str_replace(names(d), fixed("_"), ""))
            cbind(melt(d, id = "date", var = "ageb", val = "pop"),
                only.65plus = only.65plus)}))))

    # The population numbers we have are per AGEB (Área Geoestadística
    # Básica, or Basic Geostatistical Area), which are subregions of
    # delegaciones.
    message("Loading AGEB definitions")
    d = local(
       {agebs = read_sf(mexico.city.agebs.path)
        agebs = data.table(
            ageb = factor(agebs$CVEGEO),
            # For each AGEB, choose a cell of the master grid by finding the
            # cell closest to the AGEB's centroid.
            mrow = get.knnx(k = 1,
                master.grid[, .(x_sinu, y_sinu)],
                st_coordinates(st_transform(crs = crs.satellite,
                    st_centroid(agebs$geometry))))$nn.index[,1])
        stopifnot(identical(levels(agebs$ageb), levels(d$ageb)))
        merge(d, agebs, by = "ageb", all = T)})
    stopifnot(!anyNA(d))

    d = merge(d, by = c("date", "mrow"), all = T, rbindlist(lapply(
        unique(year(d$date)),
        function(y)
           {message("Loading predictions: ", y)
            p = predict.temps(y, "pred.area")[mrow %in% d$mrow]
            p[, date := as.Date(paste0(y, "-01-01")) + (yday - 1)]
            p[, yday := NULL]
            p[date %in% d$date]})))
    stopifnot(!anyNA(d))

    message("Summarizing by delegación")
    d[, by = .(only.65plus, date, delegacion = factor(substr(ageb, 1, 5))),
        .SDcols = paste0("pred.", temp.ground.vars),
        lapply(.SD, function(temp) sum(temp * pop)/sum(pop))]}
deleg.weighted.preds = pairmemo(deleg.weighted.preds, pairmemo.dir, fst = T)

per.mrow.population = function(pop.col)
   {message("Making subgrid")
    subgrid = master.grid[(in.pred.area)]
    subgrid.ps = st_sf(subgrid[, .(mrow)],
        st_sfc(crs = crs.satellite, mapply(SIMPLIFY = F, square.xyd,
            subgrid$x_sinu, subgrid$y_sinu,
            master.grid[, median(c(
                diff(sort(unique(x_sinu))),
                diff(sort(unique(y_sinu)))))])))

    message("Getting AGEBs")
    agebs = st_transform(crs = crs.satellite, read_sf(all.agebs.path))

    message("Intersecting")
    agebs$ageb.area = st_area(agebs)
    isect = st_intersection(agebs, subgrid.ps)
    isect$pop = with(isect,
        units::drop_units(get(pop.col) * st_area(isect) / ageb.area))
    setDT(isect)
    isect = isect[, by = mrow, .(pop = sum(pop))]
    isect = rbind(isect,
        data.table(mrow = setdiff(subgrid$mrow, isect$mrow), pop = 0))
    setkey(isect, mrow)
    isect}
per.mrow.population = pairmemo(per.mrow.population, pairmemo.dir, mem = T, fst = T)

square.xyd = function(x, y, d)
  # Create a polygon object representing a square from the x- and
  # y-coordinates of the center and its side length `d`.
   st_polygon(list(matrix(ncol = 2, byrow = T, c(
       x + d/2, y + d/2,    # NE
       x - d/2, y + d/2,    # NW
       x - d/2, y - d/2,    # SW
       x + d/2, y - d/2,    # SE
       x + d/2, y + d/2)))) # NE

c(master.grid, ground, stations, stns.by.dist) %<-% get.nonsatellite.data()
mrow.sets$pred.area = master.grid[, which(in.pred.area)]
mrow.sets$test = mrow.sets$pred.area[1:50]
