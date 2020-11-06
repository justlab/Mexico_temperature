# The main entry point is `save.ground`.

# The tag "XTRA" is used to mark places where I've omitted
# support for certain rare formats of the station files. By adding
# such support, one could process a few more files and hence get
# a few more observations (any number of which may not be in the
# study area, anyway).

# Warning about pressure: it isn't clear if any of the pressures
# values have been adjusted to sea level. This should be investigated
# before using pressure values.

## -----------------------------------------------------------
## * Libraries and constants
## -----------------------------------------------------------

suppressPackageStartupMessages(
   {library(data.table)
    library(readxl)
    library(XML)
    library(DBI)
    library(RSQLite)
    library(jsonlite)
    library(stringr)
    library(readr)
    library(measurements)
    library(pbapply)
    library(sf)
    library(nngeo)})

source("common.R")
source("earthdata.R")

proportion.of.day.required = .80

# In the output, all dates signify UTC-06:00 (except for Wunderground).
target.tz = "Etc/GMT+6"
  # Yes, the sign is intepreted in the opposite fashion from usual.

spanish.month.abbrs = c("ene", "feb", "mar", "abr", "may", "jun", "jul", "ago", "sep", "oct", "nov", "dic")

## ------------------------------------------------------------
## * Subroutines
## ------------------------------------------------------------

stpath = function(...) file.path(data.root, "stations", ...)

download = function(url, to, f, ..., user = NULL, password = NULL)
  # Downloads a file from `url` and saves it to `to` if there isn't
  # already a file there, then calls `f(to, ...)`.
   {to = stpath("downloads", to)
    dir.create(dirname(to), showWarnings = F, recursive = T)
    if (!file.exists(to))
        stopifnot(0 == system2("wget", c(
            url, "-O", to,
            "--no-check-certificate",
            "--user-agent", "some-program",
            (if (!is.null(user)) c("--user", user)),
            (if (!is.null(password)) c("--password", password)))))
    f(to, ...)}

slurp.zip = function(path)
  # Extracts each file in the ZIP file at `path` as a single
  # string. (Zero-length files are ignored.) The names of the
  # resulting list are set to the inner file names.
  # https://stackoverflow.com/a/56191644
   {fs = unzip(path, list = T)
    sapply(fs$Name[fs$Length > 0], simplify = F, function(x)
        paste0(collapse = '\n', read_file(unz(path, x))))}

mlr = function(...)
    regex(..., multiline = T)

if.enough.halfhourly = function(f, vals)
   {if (sum(!is.na(vals)) >= (24*2) * proportion.of.day.required)
       f(vals, na.rm = T)
    else
       NA_real_}

daily.summary = function(d, freq, variable.name)
   {d = d[!is.na(value)]
    stopifnot(!anyNA(d))

    max.vals.per.day =
       (if (freq == "hour") 24 else
        if (freq == "half hour") 24*2 else
        if (freq == "10 min") (24*60)/10 else stop())
    vals.required = ceiling(max.vals.per.day * proportion.of.day.required)

    d = d[, by = .(date, stn), eval(bquote(
      {stopifnot(.N <= max.vals.per.day)
       if (.N >= vals.required) .(
           if (variable.name == "temp.C")
               quote(list(
                    temp.C.mean = mean(value),
                    temp.C.max = max(value),
                    temp.C.min = min(value)))
           else if (variable.name == "precipitation.mm")
                quote(list(precipitation.mm.total = sum(value)))
           else
                quote(list(value.mean = mean(value))))}))]
    if ("value.mean" %in% names(d))
        setnames(d, "value.mean", paste0(variable.name, ".mean"))
    d}

combine.dailies = function(l)
    Reduce(x = l, function(x, y)
        merge(x, y, all = T, by = c("stn", "date")))

rename.cols = function(d, originals, target)
   {matches = intersect(originals, colnames(d))
    stopifnot(length(matches) <= 1)
    if (length(matches))
       setnames(d, matches, target)}

is.subset = function(needles, haystack)
    !length(setdiff(needles, haystack))

## ------------------------------------------------------------
## * Per-network functions
## ------------------------------------------------------------

## ** SIMAT

pm(get.ground.raw.simat <- function()
  # According to an email exchange Iván Gutiérrez-Avila had with
  # somebody at SIMAT, all times are in UTC-06:00.
  # Info:
  #   Air pressure: http://web.archive.org/web/20190412142523/http://www.aire.cdmx.gob.mx/descargas/datos/excel/Presion-Atmosferica.pdf
  #   Everything else: https://web.archive.org/web/20190412141446/http://www.aire.cdmx.gob.mx/descargas/datos/excel/REDMETxls.pdf
  # There's precipitation data available, but only with weekly samples:
  #   http://148.243.232.112:8080/opendata/redda/concentracion.csv - code "PP"
   {years = year(earliest.date) : latest.year
    url.root = "http://www.aire.cdmx.gob.mx/opendata"

    message("Loading SIMAT")
    obs = rbindlist(fill = T, pblapply(years, function(dyear)
       {l = c(
            list(download(
                sprintf("%s/anuales_horarios_gz/meteorología_%d.csv.gz",
                    url.root, dyear),
                sprintf("simat/meteorología_%d.csv.gz", dyear),
                fread)),
            if (dyear < 2009) list() else list(download(
              # Air-pressure observations start on this year.
                sprintf("%s/presion/PA_%d.csv", url.root, dyear),
                sprintf("simat/PA_%d.csv", dyear),
                fread)))
        for (d in l)
            setnames(d, c("datetime", "stn", "variable", "value", "unit"))
        d = rbindlist(l)
        d[, date := as.Date(
            gsub("-", "/", datetime, fixed = T),
              # Pressure files use hyphens whereas meterology files
              # use slashes.
            format = "%d/%m/%Y")]
        d[, variable :=
            ifelse(variable == "TMP", "temp.C",
            ifelse(variable == "WSP", "wind.speed.mps",
            ifelse(variable == "RH",  "relative.humidity.percent",
            ifelse(variable == "PA",  "pressure.mmHg",
                                      "other"))))]
        combine.dailies(lapply(
            setdiff(unique(d$variable), "other"),
            function(vname)
                daily.summary(d[variable == vname], "hour", vname)))}))

    stations = download(
        sprintf("%s/catalogos/cat_estacion.csv", url.root),
        "simat/cat_estacion.csv",
        fread, skip = 1, encoding = "Latin-1")
    stations = stations[, .(
        stn = cve_estac, name = nom_estac,
        lon = longitud, lat = latitud)]
    stopifnot(!anyNA(stations))
    setkey(stations, stn)

    punl(stations, obs)})

## ** UNAM

pm(get.ground.raw.unam <- function()
   {url.root = "http://www.ruoa.unam.mx/pembu/Estaciones"

    query = as.data.table(expand.grid(
        stn = c("CCA", paste0("ENP", 1:9), paste0("CCH", c("A", "N", "O", "S", "V"))),
        dyear = year(earliest.date) : latest.year))
    query = query[!(stn == "CCA" & dyear < 2007)]
      # CCA observations start on this year.

    stations = data.table()
    message("Loading UNAM")
    obs = unique(rbindlist(unlist(rec = F, pblapply(1 : nrow(query), function(qi)
       {stn.nominal = query[qi, stn]
        dyear = query[qi, dyear]
          # The station is "nominal" in the sense that files for one
          # station can have observations that are stated inside the file
          # to come from a different station.

        files = download(
            sprintf("%s/%s/datos/%s/%d_%s.zip",
                url.root, tolower(stn.nominal), stn.nominal, dyear, stn.nominal),
            sprintf("unam/%d_%s.zip", dyear, stn.nominal),
            function(path)
               {filetype = system2("file", c("-b", path), stdout = T)
                if (!grepl("Zip", filetype))
                  # We actually got an error page, but HTTP 200
                  # was returned anyway.
                    return(NULL)
                fnames = unzip(path, list = T)$Name
                if (length(fnames) > 3 && length(fnames) <= 14)
                    slurp.zip(path)
                else
                  # XTRA
                    NULL})
        if (is.null(files))
            return(list(data.table()))

        # Remove any Fortran programs that were accidentally included
        # in the archive.
        for (fname in str_subset(names(files), "\\.for\\z"))
            files[[fname]] = NULL

        # Sometimes, files seem to consist of multiple CSV
        # files squished together, so break them apart.
        files = unlist(lapply(unname(files), function(text)
            str_split(text, "\nPrograma de Estaciones")[[1]]))

        lapply(files, function(text)
           {stopifnot(str_match(text, mlr("^Horario (\\S+)"))[,2]
                == "UTC-6h")

            stn.name = str_replace_all(
                str_match(text, mlr("^Estacion (.+?)-UNAM,Ciudad de Mexico"))[,2],
                "\\s+", "")
            if (stn.name == "CCH0")
                stn.name = "CCHO"
            if (stn.name %in% c("ENP4", "ENP9"))
              # These stations sometimes have more than one
              # observation for the same time with different values
              # (e.g., compare
              # "2007/10/31 08:00:00" in
              #     2007_ENP4.zip/2007-ENP4-L1/2007-10-ENP4-L1.CSV
              # to the same row in
              #     2007_CCA.zip/2007-CCA-L1/2007-10-CCA-L1.CSV
              # .) Let's ignore them.
                return(data.table())
            latlon = str_match(text, mlr("^Lat ([0-9.]+) N, Lon ([0-9.]+) W"))
            stations <<- rbind(stations, data.table(
                stn = stn.name,
                lon = -as.numeric(latlon[,3]),
                lat = as.numeric(latlon[,2])))

            d = fread(
                str_sub(text, str_locate(text, mlr("^[0-9]{4}/"))[,1]),
                col.names = str_split(str_extract(text, mlr("^Fecha_hora,.+")), ",")[[1]],
                na.strings = "null")

            # Read dates, and drop illegal dates like September 31st.
            d[, date := as.Date(Fecha_hora)]
            d = d[!is.na(date)]

            d = d[, by = date, .(
                stn = stn.name,
                temp.C.mean = if.enough.halfhourly(mean, Temp),
                temp.C.max = if.enough.halfhourly(max, Temp),
                temp.C.min = if.enough.halfhourly(min, Temp),
                wind.speed.mps.mean = if.enough.halfhourly(mean, Rapidez_v_sostenido),
                relative.humidity.percent.mean = if.enough.halfhourly(mean, Hum_Rel),
                pressure.hPa.mean = if.enough.halfhourly(mean, Presion_bar))]
            # Drop rows that are missing on all the data columns.
            vnames = setdiff(colnames(d), c("date", "stn"))
            d[rowSums(!is.na(d[, mget(vnames)])) > 0]})}))))

    # When a single station name is associated with more than one
    # lon-lat pair, use the most common one.
    stations = (stations
        [, by = .(stn, lon, lat), .N]
        [, keyby = stn, .SD[which.max(.N)]])
    # Correct the coordinates for CCHS (CCH-SUR) by hand. These
    # numbers are from
    # https://www.ruoa.unam.mx/pembu/index.php?page=cchs
    # , but what appear to be decimal degrees there are actually
    # minutes and seconds.
    stations[.("CCHS"), c("lon", "lat") := .(
        -(99 + 11/60 + 56/(60*60)),
        19 + 18/60 + 44/(60*60))]

    punl(stations, obs)})

## ** SMN

## *** Observatories

pm(get.ground.raw.smn.observatories <- function()
  # We use the hourly file. This network also provides daily and
  # 15-minute files. However, judging from the ELEMENT-CODE
  # column of these files, they may not have the same variables
  # available.
   {message("Loading SMN observatories")

    na.value = -99999

    # Read in the giant CSV of hourly observations.
    whole = fread(cmd = paste("unzip -p", shQuote(stpath(
        "smn-raw", "OBS_2018", "Observatorios_Horarios.zip"))))
    setnames(whole, colnames(whole),
        gsub("-", "_", colnames(whole), fixed = T))
    setnames(whole,
        c("Station_ID", "ELEMENT_CODE", "YEAR_MONTH_DAY"),
        c("stn", "variable", "date"))
    whole = whole[date >= earliest.date]

    # Read dates, and drop illegal dates like September 31st.
    whole[, date := as.Date(date)]
    whole = whole[!is.na(date)]

    # Read the file of variable codes and clean up the resulting
    # variable names.
    variable.codes = fread(paste(collapse = "\n", sub("\t\t", "\t", local(
       {o = file(encoding = "Windows-1252", stpath(
            "smn-raw", "OBS_2018",
            "ELEMENTOS_HORARIOS_OBSERVATORIOS.txt"))
        x = readLines(o, warn = F)
        close(o)
        x}))))
    setkey(variable.codes, ELEMENTO)
    whole[, variable := factor(variable.codes[.(whole$variable), factor(
        ifelse(NOMBRE == "TEMPERATURA BULBO SECO °C",  "temp.C",
        ifelse(NOMBRE == "PRECIPÍTACION HORARIA mm",   "precipitation.mm",
        ifelse(NOMBRE == "HUMEDAD RELATIVA %",         "relative.humidity.percent",
        ifelse(NOMBRE == "PRESION DE LA ESTACION HPa", "pressure.hPa",
                                                       "other")))))])]

    # Collect each variable.
    obs = combine.dailies(lapply(
        setdiff(unique(whole$variable), "other"),
        function(vname)
           {d = melt(
                whole[variable == vname,
                    mget(c("stn", "date",
                        grep("VALUE", colnames(whole), val = T)))],
                id.vars = c("date", "stn"),
                variable.name = "time",
                value.name = "value")
            d[value == na.value, value := NA]
            daily.summary(d, "hour", vname)}))

    # Collect stations.
    stations = as.data.table(read_excel(stpath(
        "smn-raw", "OBS_2018", "Claves_Observatorios.xlsx")))
    stations = stations[!is.na(CLAVE), .(
        stn = as.integer(CLAVE),
        lon = -(LONG + LONM/60 + ifelse(is.na(LONS), 0, LONS/(60*60))),
        lat = LATG + LATM/60 + ifelse(is.na(LATS), 0, LATS/(60*60)))]
    
    punl(stations, obs)})

get.ground.raw.smn.observatories.new <- function()
  # Incomplete code to collect newer versions of the SMN observatories
  # files from an FTP site.
   {dl = function(fname, f, ...)
        download(
            paste0("ftp://200.4.8.36:/BD-Climatologica/Observatorios/",
               fname),
            file.path("smn-observatories", fname),
            user = "ftpsmn.conagua",
            password = Sys.getenv("JUSTLAB_MEXICO_TEMPERATURE_SMN_FTP_PASSWORD"),
            f = f, ...)

    dl("ELEMENTOS_HORARIOS_OBSERVATORIOS.txt", message)
    dl("Claves_Observatorios.xlsx", message)
    dl("OBS_HLY.zip", message)}

## *** ESIMEs and EMAs

pm(get.ground.raw.smn.esimes <- function()
    es.stations(process.es.observations(read.es())))

pm(get.ground.raw.smn.emas <- function()
    es.stations(emas = T,
        process.es.observations(emas = T, n.jobs = 8, read.es(emas = T))))

read.es = function(emas = F)
   {message("Loading SMN ", (if (emas) "EMAs" else "ESIMEs"))

    fnames = list.files(
       (if (emas)
           stpath("smn-emas-csv") else
           stpath("smn-raw", "ESIMEs_2018")),
       full.names = T)

    if (!emas)
       {fnames = grep(".xls", fnames, val = T)
          # XTRA: CSV ESIMEs files
        fnames = grep("ESIME Guanajuato 2009", invert = T, fnames, val = T)
          # XTRA: this file leads to an error due to a row that's
          # empty except for the first cell, which is "faltantes"
        fnames = grep("ESIME Valladolid 2009", invert = T, fnames, val = T)}
          # XTRA: the last row is empty except for the first
          # cell, which is "Faltantes de oct, nov y dic"

    ds = pblapply(fnames, function(fname)
       {if (endsWith(fname, ".csv"))
           {d = fread(fname, na.strings = "")
            if (nrow(d) && is.na(d[1, 1]))
              # We mistook a row of units for a row of values.
              # Reread the file, adding the units to the column names.
               {units = sapply(d[1,], as.character)
                d = as.data.table(fread(fname, na.strings = "",
                    col.names = ifelse(is.na(units),
                        colnames(d),
                        sprintf("%s(%s)", colnames(d), units)),
                    skip = 2))}}
        else
           {d = tryCatch(read_excel(fname),
                error = function(e)
                    list(fname = fname, em = conditionMessage(e)))
            if ("em" %in% names(d))
                return(d)
            setDT(d)
            if (nrow(d) && is.na(d[1, 1]))
              # We mistook a row of units for a row of values.
              # Reread the file, adding the units to the column names
              # to be more like the other files.
               {units = sapply(d[1,], as.character)
                d = as.data.table(read_excel(fname,
                    col_names = ifelse(is.na(units),
                        colnames(d),
                        sprintf("%s(%s)", colnames(d), units)),
                    range = cell_limits(c(3, 1), c(NA, length(colnames(d))))))}}
        d})

    names(ds) = sapply(fnames, basename)
    ds}

process.es.observations = function(ds, emas = F, n.jobs = NULL)
   {di = 0
    d = rbindlist(fill = T, Filter(nrow, pblapply(ds, cl = n.jobs, function(d)
       {di <<- di + 1
        #message(di)
        fname = names(ds)[di]

        # Skip empty files.
        if (!nrow(d))
            return(d)

        if (fname %in% c(
                "ESIME Aguascalientes 2006.xlsx",
                  # This file has repeated DateTimes with
                  # distinct temperatures.
                "ESIME Tepehuanes 2017.xlsx"))
                  # This file only has measurements from earlier
                  # years, which don't match what's in the
                  # earlier years' files.
            return(data.table())

        if (!emas)
          # For ESIMEs, skip any file that appears to contains more
          # than one station, because we rely on the file name to
          # determine the station name.
           {stn.col = str_subset(colnames(d), "^(Estaci.+?n|Station)$")
            stopifnot(length(stn.col) == 1)
            if (length(unique(d[[stn.col]])) > 1)
                return(data.table())}

        if (emas &&
                !is.subset(c("fecha", "TempAire", "RapViento", "HumRelativa", "PresBarometric", "Precipitacion", "nombre_estacion"), colnames(d)) &&
                !is.subset(c("Date", "Time", "AvgTemp(C)", "AvgBP(mbar)", "AvgRh(%)", "WSK(kph)", "WSMK(kph)", "Rain(mm)"), colnames(d)))
         # XTRA: other EMAs formats aren't considered.
            return(data.table())

        d = copy(d)

        # Set the DateTime column.
        rename.cols(d, c("Fecha-Tiempo", "fecha"), "DateTime")
        if ("Date" %in% colnames(d) & "Time" %in% colnames(d))
           {dparts = str_match(tolower(d$Date), "(\\d{4}) ([a-z]{3}) (\\d\\d?)")
            d[, DateTime := as.POSIXct(
                tz = "UTC",
                format = "%Y %m %d %I:%M:%S %p",
                paste(
                    dparts[,2],
                    match(dparts[,3], spanish.month.abbrs),
                    dparts[,4],
                    d$Time))]}
        stopifnot("DateTime" %in% colnames(d))
        d = d[!is.na(DateTime)]
        # Skip any file that has become empty.
        if (!nrow(d))
            return(d)

        # Normalize station names.
        rename.cols(d, "nombre_estacion", "stn")
        if (emas)
           {if (!("stn" %in% colnames(d)))
                d[, stn := str_replace(tolower(fname),
                   " +((emz|mrab|mzab|abjn|myjn|jlag) *)?([a-z][a-z]? *)?[0-9][0-9][an]? *\\.csv",
                   "")]}
        else
            d[, stn := str_match(fname, "ESIME (.+) \\d{4} *\\.")[,2]]

        # Skip an EMAs station with many repeated datetimes.
        if (d$stn[1] == "acapo")
            return(data.table())

        ifcol = function(...)
            for (col in c(...))
                if (col %in% colnames(d))
                   {v = d[[col]]
                    if (is.character(v))
                        v = as.numeric(str_replace(v, fixed(","), "."))
                    return(v)}
        d = d[, .(
            stn,
            datetime =
              # Iván Gutiérrez-Avila says all timestamps are in UTC.
              # This is plausible given the observed temperatures.
               {if (is.character(DateTime))
                  {v = str_replace_all(tolower(DateTime), fixed("."), "")
                   as.POSIXct(v,
                       format = (if (any(str_detect(v, "\\d{8}")))
                           "%Y%m%d %H%M" else
                           paste(
                               (if (any(str_detect(v, "/(?:[2-3][0-9]|1[3-9])/")))
                                   "%m/%d/%Y" else
                                   "%d/%m/%Y"),
                               (if (any(str_detect(v, " pm")))
                                   "%I:%M:%S %p" else
                                   "%H:%M:%S"))),
                       tz = "UTC")}
                else
                   {stopifnot(attr(DateTime, "tz") == "UTC")
                    DateTime}},
            temp.C = ifcol(
                "ATC(C)", "TempAire", "AvgTemp(C)"),
            relative.humidity.percent = ifcol(
                "RH %", "Rh(%)", "HumRelativa", "AvgRh(%)"),
            precipitation.mm = ifcol(
                "Rain(mm) 10Min", "Rain(mm)", "Precipitacion"),
              # XTRA: some other ESIMEs column names that are infrequent:
              # "Rain(mm) 1Hr", "Rain(mm) 24hrs", "Rain(mm) 3Hrs", "Rain(mm) 6Hrs"
              # The unit for "Precipitacion" is a guess based on other
              # EMAs files.
            wind.speed.mps = ifcol("WS(m/s)"),
            wind.speed.kmph = ifcol("WSK(kph)", "RapViento"),
              # The unit for RapViento is a guess based on the units
              # for other EMAs files.
              # XTRA: there are also some files with wind speed *u* vs.
              # wind speed *v* ("AvgWSU(m/s)", "AvgWSV(m/s)").
            pressure.hPa = ifcol("BP(mbar)", "AvgBP(mbar)", "PresBarometric"))]

        # Keep only observations at the 10-minute intervals. I
        # see other observations in only a few cases, and they're
        # at weird times like 5 seconds past a 10-minute
        # interval.
        d[minute(datetime) %% 10 == 0 & second(datetime) == 0]})))

    message("Intermediate processing")
    d = unique(d)
    # Drop rows that are missing on all the data columns.
    vnames = setdiff(colnames(d), c("date", "datetime", "stn"))
    d = d[rowSums(!is.na(d[, mget(vnames)])) > 0]
    # Drop a single pesky partly duplicated observation.
    d = d[!(stn == "la rumor" &
       datetime == as.POSIXct(tz = "UTC", "2006-07-21 22:00:00") &
       is.na(relative.humidity.percent))]
    # We should now have only one observation per station and time.
    stopifnot(all(d[, .N, by = .(stn, datetime)]$N == 1))
    # Add dates.
    d[, date := as.Date(datetime, tz = target.tz)]

    message("Summarizing")
    combine.dailies(lapply(
        setdiff(colnames(d), c("stn", "date", "datetime")),
        function(vname)
            daily.summary(
                d[, .(stn, date, value = get(vname))],
                "10 min",
                  # XTRA: some files have some variables present
                  # only at a coarser interval, such as 1 hour.
                vname)))}

es.stations = function(obs, emas = F)
   {stations = as.data.table(download(
        "http://web.archive.org/web/20180809001831id_/http://smn1.conagua.gob.mx/emas/catalogoa.html",
        "smn_es_stations.html",
        function(fname)
            readHTMLTable(fname,
                encoding = "Windows-1258",
                which = (if (emas) 5 else 6),
                header = T, stringsAsFactors = F)))

    coord = function(s)
       {s = str_replace_all(s, fixed(","), ".")
        e = apply(MARGIN = 2, FUN = as.numeric,
            str_match(s, "(\\d+)[°º]\\s*(\\d+)['´]\\s*([0-9.]+)")[, -1])
        e[,1] + e[,2]/60 + e[,3]/(60*60)}
    stations = stations[, .(
        stn = NOMBRE,
        lon = -coord(`Longitud (Oeste)`),
        lat = coord(`Latitud (Norte)`))]
    # When stations are duplicated, keep only the first.
    stations = stations[, by = stn, head(.SD, 1)]
    stopifnot(!anyNA(stations))

    l = list(seen = unique(obs$stn), known = stations$stn)
    l$seen.orig = l$seen

    for (x in c("seen", "known"))
        {l[[x]] = str_replace(l[[x]], "\\xd1", "")
         l[[x]] = iconv(tolower(l[[x]]), to = "ASCII//TRANSLIT")

         replacements = c(
             "\\s+" = " ",
             "^c\\." = "canon ",
             "\\bcanon del?\\b" = "canon",
             "\\bbh\\b" = "bahia",
             "\\bbahia del?( los)?\\b" = "bahia",
             "\\bang[a-z]{0,2}$" = "angeles",
             "\\bs luc(as)?\\b" = "san lucas",
             "^cer " = "cerro ",
             "^nev " = "nevado ",
             "^pq " = "parque ",
             "^psa " = "presa ",
             "^r " = "rio ",
             "^s " = "san ",
             "^v " = "villa ",
             "\\bnevado de " = "nevado ",
             ",.+" = "",
             "\\s+\\(.+?\\)" = "",
             "\\b(ciudad\\b|cd\\b\\.?)" = "",
             "\\bgpe$" = "",
             " ipn$" = "",
             "\\bencb" = "escuela nacional de ciencias biologicas",
             "\\bescuela nacional de ciencias biol.?gicas\\z" = "escuela nacional de ciencias biologicas i",
             "2$" = "ii",
             "^utt$" = "universidad tecnologica de tecamachalco",
             "^gust.*" = "gustavo diaz ordaz",
             "^jo ma mor.*" = "jose maria morelos",
             "\\bnacional\\b" = "",
             "^cuatroc$" = "cuatro cienegas",
             "\\bp\\.l.?pez zamora" = "presa emilio lopez zamora",
             ".+?emil.+" = "presa emilio lopez zamora",
             ".+?abelar.+" = "presa abelardo l. rodriguez",
             "\\bcumbres de mty - el diente" = "el diente",
             "\\bvilla carranza" = "venustiano carranza",
             "\\bnochix\\b" = "nochistlan",
             "^presa cuchi" = "presa el cuchillo",
             "(presa )?\\bla cang.*" = "presa cangrejera",
             "\\bjuani\\S*" = "juanico",
             "\\btantak\\b" = "tantaquin",
             "\\btlapa\\b.*" = "tlapa de comonfort",
             "\\bpantanos de centla\\b" = "centla")
         for (r in names(replacements))
             l[[x]] = str_replace_all(l[[x]], r, replacements[[r]])

         l[[x]] = trimws(str_replace_all(l[[x]], "\\s+", " "))}

    stopifnot(!anyDuplicated(l$known))

    dists = adist(l$seen, l$known, ignore.case = T)
    match.i = apply(dists, 1, which.min)
    # If a seen station is an exact prefix of exactly one known
    # station, that takes precedence over any match with high
    # Levenshtein distance.
    switch.threshold = 1
    match.i = sapply(seq_along(l$seen), function(i)
       {if (min(dists[i,]) > switch.threshold)
           {prefix.of = startsWith(l$known, l$seen[i])
            if (sum(prefix.of) == 1)
                return(which(prefix.of))}
        match.i[i]})
    # Manually unmatch stations that, in my eyes, don't have
    # exactly one good match.
    no.matches = c(
        "^cca$", "^leon$", "\\bvallarta\\b", "iztacihuatl",
        "^matat$", "\\bloma\\b", "^montemorelos$",
        "^tampico madero$", "^queretaro conagua$", "^merida conagua$")
    for (x in no.matches)
        match.i[str_detect(l$seen, x)] = NA

    # Print matches, for debugging.
    #print(rbind(
    #    unique(obs$stn),
    #    tolower(stations[match.i, stn]),
    #    l$seen,
    #    l$known[match.i]))

    obs.matches = match.i[match(obs$stn, l$seen.orig)]
    # Drop observations with no matching station.
    obs = obs[!is.na(obs.matches)]
    obs[, stn := stations[na.omit(obs.matches), stn]]

    punl(stations, obs)}

## ** Wunderground

# Slippage: Wunderground's dates are in America/Mexico_City,
# not `target.tz`.
#
# Further slippage: what I interpret as the mean pressure may not
# actually be the mean.

get.ground.raw.wunderground = function()
   {db = dbConnect(SQLite(),
        stpath("wunderground-daily-mexico.sqlite"))
    on.exit(dbDisconnect(db))
    stations = as.data.table(dbGetQuery(db, "select
            stn,
            lon,
            lat,
            elevation as 'stn.elevation',
            height,
            surface_type as surface,
            state,
            city,
            neighborhood,
            station_type as 'station.type',
            software
        from Stations"))
    obs = as.data.table(rbind(
        dbGetQuery(db, sprintf("select %s from Daily_Old_API",
           "stn,
            date,
            max_temperature as 'temp.F.max',
            min_temperature as 'temp.F.min',
            temperature as 'temp.F.mean',
            precip_today as `precipitation.inch.total`,
            wind_speed as `wind.speed.miles.per.hour.mean`,
            pressure as `pressure.inHg.mean`")),
        dbGetQuery(db, sprintf("select %s from Daily",
           "stn,
            date,
            imperial_tempHigh as 'temp.F.max',
            imperial_tempLow as 'temp.F.min',
            imperial_tempAvg as 'temp.F.mean',
            imperial_precipTotal as `precipitation.inch.total`,
            imperial_windspeedAvg as `wind.speed.miles.per.hour.mean`,
            imperial_pressureTrend as `pressure.inHg.mean`"))))
    obs[, date := as.Date(as.character(date), format = "%Y%m%d")]
    stopifnot(!anyNA(obs$date))
    punl(stations, obs)}

## ------------------------------------------------------------
## * All-network functions
## ------------------------------------------------------------

save.ground = function()
   {on.exit(close(con))
    con = gzfile(ground.json.path)
    cat(file = con, toJSON(
        dataframe = "columns", digits = NA, na = "null",
        do.call(filter.raw, get.ground.raw())))}

get.ground.raw = function()
   {networks = list(
        simat = get.ground.raw.simat(),
        unam = get.ground.raw.unam(),
        smno = get.ground.raw.smn.observatories(),
        esimes = get.ground.raw.smn.esimes(),
        emas = get.ground.raw.smn.emas(),
        wunderground = get.ground.raw.wunderground())

    # Add network identifiers.
    for (net in names(networks))
       {stations = copy(networks[[net]]$stations)
        obs = copy(networks[[net]]$obs)
        stations[, network := net]
        for (x in list(stations, obs))
            x[, stn := paste(stn, net)]
        networks[[net]] = punl(stations, obs)}

    # Combine all networks.
    stations = rbindlist(fill = T, lapply(networks, function(x)
        cbind(x$stations[, .(stn, network, lon, lat)],
            if (x$stations$network[1] == "wunderground")
              # Include all the other columns, too.
                  x$stations[, .SD, .SDcols =
                      !c("stn", "network", "lon", "lat")])))
    obs = rbindlist(fill = T, lapply(networks, function(x)
        x$obs))

    # Drop stations with no observations.
    stations = stations[stn %in% obs$stn]

    stopifnot(!anyNA(obs$date))
    stopifnot(!anyNA(obs$stn))
    stopifnot(is.subset(unique(obs$stn), stations$stn))
    stopifnot(nrow(unique(obs[, .(stn, date)])) == nrow(obs))

    # Convert units.
    for (unit in list(
            list(from = "F", to = "C"),
            list(from = "mmHg", to = "hPa"),
            list(from = "inHg", to = "hPa"),
            list(from = "kmph", to = "mps"),
            list(from = "miles.per.hour", to = "mps"),
            list(from = "inch", to = "mm")))
       {reg = sprintf("\\.%s\\.(mean|min|max|total)$", unit$from)
        cn = copy(colnames(obs))
        for (col.from in cn)
            if (str_detect(col.from, reg))
               {col.to = str_replace(col.from, reg,
                    sprintf(".%s.\\1", unit$to))
                obs[, (col.to) := ifelse(is.na(get(col.to)),
                    (if (unit$from == "kmph" && unit$to == "mps")
                            conv_multiunit(get(col.from), "km / hr", "m / sec") else
                        if (unit$from == "miles.per.hour" && unit$to == "mps")
                            conv_multiunit(get(col.from), "mi / hr", "m / sec") else
                            conv_unit(get(col.from), unit$from, unit$to)),
                    get(col.to))]
                obs[, (col.from) := NULL]}}

    punl(stations, obs)}

min.obs = 100
possible.duplicate.dist.meters = 2000
min.common.days = 30
max.proportion.equal = 1/10

deviance.quantile = .95
idw.maxdist.meters = 30e3
idw.elev.thresh.meters = 500

# https://en.wikipedia.org/wiki/Climate_of_Mexico#Weather_records
extreme.hi.temp.C = 53
extreme.lo.temp.C = -30
extreme.precipitation.mm = 1634
# https://en.wikipedia.org/wiki/Wind_speed
extreme.wind.speed.mps = 114

filter.raw = function(stations, obs)
   {status = function()
       message(sprintf("Now at %s observations from %s stations",
           format(nrow(obs), big.mark = ","),
           format(nrow(stations), big.mark = ",")))
    status()

    message("Narrowing to study area")
    stations = stations[in.study.area(lon, lat)]
    obs = obs[stn %in% stations$stn]
    status()

    message("Removing stations with only a few observations")
    obs = obs[, by = stn, if (.N >= min.obs) .SD]
    stations = stations[stn %in% unique(obs$stn)]
    status()

    message("Removing station-days with absurd values")
    local(
       {absurd = obs[, .(
            hot.max = temp.C.max >= extreme.hi.temp.C,
            hot.avg = temp.C.mean >= extreme.hi.temp.C,
            hot.min = temp.C.min >= extreme.hi.temp.C,
            cold.max = temp.C.max <= extreme.lo.temp.C,
            cold.avg = temp.C.mean <= extreme.lo.temp.C,
            cold.min = temp.C.min <= extreme.lo.temp.C,
            inconsistent = (temp.C.max < temp.C.min  |
                temp.C.mean < temp.C.min |
                temp.C.mean > temp.C.max),
            wet = precipitation.mm.total >= extreme.precipitation.mm,
            dry = precipitation.mm.total < 0,
            fast = wind.speed.mps.mean >= extreme.wind.speed.mps,
            slow = wind.speed.mps.mean < 0)]
        absurd = absurd[, lapply(.SD, function(v)
            ifelse(is.na(v), F, v))]
        print(colSums(absurd))
        obs <<- obs[rowSums(absurd) == 0]
        stations <<- stations[stn %in% unique(obs$stn)]})
    status()

    # Compute per-station temperature precisions.
    local(
       {close.to = function(values, comparison)
            mean(abs(values - comparison)) < .0001
        precs = obs[, keyby = stn, .(prec =
           {C = na.omit(unique(c(temp.C.max, temp.C.min)))
            f = conv_unit(C, "C", "F")
            if (close.to(C, round(C)))
                "1 C"
            else if (close.to(f, round(f)))
                "1 F"
            else if (close.to(C, round(C, 1)))
                "0.1 C"
            else if (close.to(f, round(f, 1)))
                "0.1 F"
            else if (close.to(C, round(C, 2)))
                "0.01 C"
            else if (close.to(f, round(f, 2)))
                "0.01 F"
            else
                "other"})]
        stations[, temp.prec := precs[.(stations$stn), prec]]})

    # Get the station distance matrix.
    stdist = units::drop_units(st_distance(st_as_sf(
       stations,
       coords = c("lon", "lat"), crs = crs.lonlat)))

    message("Checking for duplicate stations (by max temperature)")
    groups = cutree(h = possible.duplicate.dist.meters,
        hclust(as.dist(stdist), method = "complete"))
    for (gn in unique(groups))
       {stns = stations[groups == gn, stn]
        if (length(stns) == 1)
            next

        for (pair in combn(stns, 2, simp = F))
           {common.dates = obs[by = date,
                stn %in% pair & !is.na(temp.C.max) & !is.na(temp.C.min),
                if (.N == length(pair)) 1]$date
            if (length(common.dates) < min.common.days)
                next

            # Collect max temperatures.
            temps = obs[by = date,
                stn %in% pair & date %in% common.dates,
                .(temp.max = temp.C.max)]
            # Round them to the precision of the less precise of
            # the two stations.
            precs = stations[stn %in% pair, temp.prec]
            temps[, temp.max := (
                if (any(precs == "0.1 C"))
                    round(10 * temp.max)
                else if (any(precs == "0.1 F"))
                    round(10 * conv_unit(temp.max, "C", "F"))
                else if (all(precs == "0.01 C"))
                    round(100 * temp.max)
                else stop())]

            # Are the rounded temperatures often equal?
            proportion.equal = (temps
                [, by = date, length(unique(temp.max)) == 1]
                [, mean(V1)])
            stopifnot(proportion.equal < max.proportion.equal)}}

    message("Getting station elevations")
    stations[, elev := extract(get.elevation(),
        st_as_sf(.SD, coords = c("lon", "lat"), crs = crs.lonlat))]
    setkey(stations, stn)
    obs[, elev := stations[.(obs$stn), elev]]

    message("Removing deviant stations (by temperature)")
    # Compute the squared differences between observed temperatures
    # and IDW interpolations from other stations. Throw away the
    # stations with the greatest such differences.
    stns.deviant = character()
    local(
       {idw.tables = repeated.idw.tables(
            locations = st_coordinates(st_transform(crs = crs.mexico.city,
                st_as_sf(crs = crs.lonlat,
                    stations, coords = c("lon", "lat")))),
            maxdist = idw.maxdist.meters,
            source.subsetter = function(i)
                # Don't interpolate a station with itself.
                (1 : nrow(stations)) != i &
                # Only use stations within an elevation threshold.
                abs(stations[i, elev] - stations$elev)
                    < idw.elev.thresh.meters)
        for (temp.var in paste0("temp.C.", c("max", "min", "mean")))
           {obs[, idw := repeated.idw(
                tables = idw.tables,
                li = match(stn, stations$stn),
                group = date,
                outcome = get(temp.var),
                progress = T)]
            diffs = obs[!is.na(get(temp.var)) & !is.na(idw), by = stn,
               {sqd = (get(temp.var) - idw)^2
                .(sqd.mean = mean(sqd), sqd.max = max(sqd))}]
            stns.deviant <<- unique(c(stns.deviant,
                diffs[j = stn,
                    sqd.mean > quantile(sqd.mean, deviance.quantile) |
                    sqd.max > quantile(sqd.max, deviance.quantile)]))}})
    obs[, idw := NULL]
    stations = stations[!(stn %in% stns.deviant)]
    obs = obs[!(stn %in% stns.deviant)]
    status()

    # Replace station identifiers with simple integers, and put
    # stations and observations in a logical order.
    stations[, name := stn]
    stations = stations[order(network, name, lon, lat)]
    stations[, stn := .I]
    setkey(stations, stn)
    obs[, stn := stations[, keyby = name, stn][.(obs$stn), stn]]
    obs = obs[order(stn, date)]

    punl(stations, obs)}
