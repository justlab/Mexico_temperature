suppressPackageStartupMessages(
   {library(data.table)
    library(FNN)
    library(lme4)
    library(ape)
    library(httr)
    library(sf)
    library(future.apply)
    library(zeallot)
    library(caret)
    library(stringr)
    library(raster)
    library(hdf5r)})

source("common.R")

c(full.satellite.grid, get.satellite.data, get.elevation) %<-% local(
   {source("earthdata.R")
    list(full.satellite.grid, get.satellite.data, get.elevation)})

all.agebs.year = 2010L

plan(multiprocess)

n.folds = 10

temp.ground.vars = c(
    "ground.temp.lo", "ground.temp.mean", "ground.temp.hi")
nontemp.ground.vars = c(
    "wind.speed.mean")

available.years = year(earliest.date) : latest.year

pm(mem = T,
get.nonsatellite.data <- function()
   {message("Loading master grid")
    master.grid <<- full.satellite.grid()[in.study.area(lon, lat)]
    stopifnot(nrow(unique(master.grid[, .(lon, lat)])) == nrow(master.grid))
    master.grid[, x_sinu_f := floor(x_sinu)]
    master.grid[, y_sinu_f := floor(y_sinu)]
    setkey(master.grid, x_sinu_f, y_sinu_f)
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

    local(
       {message("Loading elevation")
        elev = get.elevation()
        message("Extracting")
        mg = st_as_sf(master.grid, coords = c("lon", "lat"), crs = crs.lonlat)
        master.grid[, elevation := extract(elev, mg)]})

    message("Saving master grid as HDF5")
    h5 = H5File$new(file.path(data.root, "master_grid.h5"), mode = "w")
    h5[["master_grid"]] = transform(as.data.frame(master.grid), region =
        factor(ifelse(is.na(region), "none", as.character(region))))
          # hdf5r can't handle factors with NAs. https://github.com/hhoeflin/hdf5r/issues/29
    h5$close_all()

    message("Loading data from ground stations")
    stations = copy(get.ground()$stations)
    ground = copy(get.ground()$obs)

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

    message("Master grid and ground data loaded")
    list(master.grid, ground, stations, stns.by.dist)})

set.mrows = function(d, longitude.col, latitude.col)
  # Adds to the given data table a column `mrow` that specifies
  # the corresponding row of `master.grid`.
   {x = d[[longitude.col]]
    y = d[[latitude.col]]
    stopifnot(all(in.study.area(x, y)))
    set(d, j = "mrow", value =
        get.knnx(master.grid[, .(lon, lat)], cbind(x, y), k = 1)$nn.index[,1])}

mrow.sets = list()

pm(mem = T, fst = T,
model.dataset <- function(the.year, mrow.set = NULL, nonmissing.ground.temp = F)
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
        wunder = stations[.(d$stn), network] == "wunderground",
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
    d})

impute.nontemp.ground.vars = function(d.orig, fold.i, progress = F, train.wunder = T)
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
                .(yday, v = get(vname), fold, wunder)]
            other.vs[[the.stn]] = rep(NA_real_, max(d$yday))
            if (train.wunder || !piece$wunder[1])
              # When we don't want to train on Wunderground stations,
              # keep the NAs so none of their values will be used.
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

pm(mem = T, fst = T,
run.cv <- function(the.year, dvname, train.wunder = T)
  # Under cross-validation, predict ground temperature using
  # satellite temperature over the given year.
   {d.master = copy(model.dataset(the.year, nonmissing.ground.temp = T))
    setnames(d.master, dvname, "ground.temp")
    d.master[, setdiff(temp.ground.vars, dvname) := NULL]
    message("run.cv: ", the.year, " ", dvname)

    bar = txtProgressBar(min = 0, max = n.folds, style = 3)
    for (fold.i in 1 : n.folds)
       {d = impute.nontemp.ground.vars(d.master, fold.i,
            train.wunder = train.wunder)
        f.pred = train.model(d[fold != fold.i &
            (train.wunder | !wunder)])
        d.master[fold == fold.i, pred := f.pred(d[fold == fold.i])]
        setTxtProgressBar(bar, fold.i)}
    close(bar)

    cbind(d.master, year = the.year, dv = dvname)})

multi.run.cv = function(years, train.wunder = T)
  # Run cross-validation for each outcome in each of the given years,
  # and combine all the results into one big data.table.
  {args = expand.grid(the.year = years, dv = temp.ground.vars,
       stringsAsFactors = F)
   rbindlist(lapply(1 : nrow(args), function(i)
       run.cv(args[i, "the.year"], args[i, "dv"],
           train.wunder = train.wunder)))}

summarize.cv.results = function(multirun.output, test.wunder = T)
   {d = multirun.output[fold != -1 & (test.wunder | !wunder)]
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

if (!exists("predict.temps.cache")) predict.temps.cache = list()

predict.temps = function(the.year, mrow.set)
  # Predict a low, mean, and high temperature for every cell in
  # `mrow.set` and day in `the.year`.
   {if (mrow.set == "pred.area")
       {pred.dir = file.path(data.root, "predictions")
        dir.create(pred.dir, showWarnings = F)
        path = file.path(pred.dir, paste0(the.year, ".h5"))
        if (path %in% names(predict.temps.cache))
            return(predict.temps.cache[[path]])}

    if (mrow.set != "pred.area" || !file.exists(path)) local(
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
            head(.SD, 1)]

        if (mrow.set != "pred.area")
            return(d)

        message("Writing HDF5")
        h5 = H5File$new(path, mode = "w")
        d = reshape2::acast(melt(d, id.vars = c("mrow", "yday")),
            mrow ~ yday ~ variable)
        h5[["data"]] = d
        h5attr(h5[["data"]], "dimensions") = c(
            "mrow", "date", "variable")
        dl = h5$create_group("dimension_labels")
        dl[["mrow"]] = as.integer(dimnames(d)[[1]])
        dl[["date"]] = as.character(
                as.Date(paste0(the.year - 1, "-12-31")) +
                as.integer(dimnames(d)[[2]]))
        dl[["variable"]] = dimnames(d)[[3]]
        h5$close_all()})

    message("Reading HDF5")
    h5 = H5File$new(path, mode = "r")
    d = h5[["data"]][,,]
    dimnames(d) = sapply(simplify = F,
        h5attr(h5[["data"]], "dimensions"),
        function(k) h5[["dimension_labels"]][[k]][])
    h5$close_all()
    d = dcast(as.data.table(melt(d)),
        mrow + date ~ variable)
    d[, date := yday(as.Date(date))]
    setnames(d, "date", "yday")
    setkey(d, mrow, yday)
    predict.temps.cache[[path]] <<- d

    d}

predict.temps.at = function(fname, date.col, lon.col, lat.col)
   {d = fread(fname)[, mget(c(date.col, lon.col, lat.col))]
    setnames(d, c("date", "lon", "lat"))
    set.mrows(d, "lon", "lat")
    d[, yday := yday(date)]
    message("Getting model predictions")
    d[, by = .(y = year(date)), paste0("model.", temp.ground.vars) :=
        predict.temps(y, "pred.area")[.(.SD$mrow, .SD$yday),
            mget(paste0("pred.", temp.ground.vars))]]

    # Also find the mean of all SIMAT stations for each day.
    simat.daily.means = with(get.ground(), obs[
        stn %in% stations[network == "simat", stn],
        keyby = date,
        mean(temp.C.mean, na.rm = T)])
    d[, simat.ground.temp.mean := simat.daily.means[.(d$date), V1]]

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
        "simat.ground.temp.mean",
        "precipitation.mm"))]}

pm(mem = T, fst = T,
per.mrow.population <- function(pop.col)
   {message("Making subgrid")
    subgrid = master.grid[(in.pred.area)]
    subgrid.ps = st_sf(subgrid[, .(mrow)],
        st_sfc(crs = crs.satellite, mapply(SIMPLIFY = F, square.xyd,
            subgrid$x_sinu, subgrid$y_sinu,
            master.grid[, median(c(
                diff(sort(unique(x_sinu))),
                diff(sort(unique(y_sinu)))))])))

    message("Getting AGEBs")
    path = file.path(data.root, "geography", paste0("agebs_", all.agebs.year))
    if (!length(list.files(path)))
       {message("Downloading")
        r = GET("http://datamx.io/dataset/a3d477e3-573a-408c-85aa-bcc6e3a8d714/resource/60aa0982-b8fa-4fb7-9e48-c3c45341c972/download/agebmexico2010.zip")
        stop_for_status(r)
        writeBin(content(r, "raw"), file.path(path, "download.zip"))
        unzip(file.path(path, "download.zip"), exdir = path)
        unlink(file.path(path, "download.zip"))
        for (item in list.files(path, full = T))
            file.rename(item, file.path(path, sprintf("agebs_%d.%s",
                all.agebs.year, tools::file_ext(item))))}
    message("Reading")
    agebs = st_transform(crs = crs.satellite, read_sf(path))

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
    isect})

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
