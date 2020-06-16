suppressPackageStartupMessages(
   {library(data.table)
    library(sf)})

source("common.R")

mexico.city.agebs.path = "~/Jdrive/PM/Just_Lab/projects/airmex/data/gis/gisdata/AGEBS_CDMX_2010.shp"
population.path.fmt = "~/Jdrive/PM/Just_Lab/projects/airmex/data/population/%s_cuadratic_csv"

pm(fst = T,
deleg.weighted.preds <- function()
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
        lapply(.SD, function(temp) sum(temp * pop)/sum(pop))]})
