# The tag "XTRA" is used to mark places where I've omitted
# support for certain rare formats of the station files. By adding
# such support, one could process a few more files and hence get
# a few more observations.

suppressPackageStartupMessages(
   {library(data.table)
    library(readxl)
    library(XML)
    library(DBI)
    library(stringr)
    library(measurements)
    library(pbapply)})

source("../Just_universal/code/pairmemo.R")
pairmemo.dir = "/data-belle/Mexico_temperature/pairmemo"
source("../Just_universal/code/punl.R")

data.root = "/data-union/Mexico_temperature/stations"
wunderground.db.path = "/data-belle/wunderground/wunderground-daily-mexico.sqlite"

earliest.date = "2003-01-01"
  # The earliest date we're interested in.
hourly.vals.required = 18  # I.e., 3/4 of the day.

# In the output, all dates signify UTC-06:00 (except for Wunderground).
target.tz = "Etc/GMT+6"
  # Yes, the sign is intepreted in the opposite fashion from usual.

## ------------------------------------------------------------
## * Subroutines
## ------------------------------------------------------------

download = function(url, to, f, ...)
  # Downloads a file from `url` and saves it to `to` if there isn't
  # already a file there, then calls `f(to, ...)`.
   {to = file.path(data.root, "downloads", to)
    dir.create(dirname(to), showWarnings = F, recursive = T)
    if (!file.exists(to))
        stopifnot(0 == system2("wget", c(
            url, "-O", to,
            "--no-check-certificate",
            "--user-agent", "some-program")))
    f(to, ...)}

slurp.zip = function(path)
  # Extracts each file in the zip file at `path` as a single
  # string. The names of the resulting list are set to the inner
  # file names.
   {lines = system2("unzip", c("-c", path), stdout = T)
    is.sep = str_detect(lines, "^ (?: inflating|extracting): ")
    chunks = lapply(
        split(lines[!is.sep], cumsum(is.sep)[!is.sep])[-1],
        function(chunk) paste(chunk, collapse = "\n"))
    fnames = str_match(lines[is.sep], "^ (?: inflating|extracting): (.+)  $")
    stopifnot(!anyNA(fnames))
    names(chunks) = fnames[,2]
    chunks}

mlr = function(...)
    regex(..., multiline = T)

if.enough.halfhourly = function(f, vals)
   {if (sum(!is.na(vals)) >= 2*hourly.vals.required)
       f(vals, na.rm = T)
    else
       NA_real_}

daily.summary = function(d, freq, variable.name)
   {d = d[!is.na(value)]
    stopifnot(!anyNA(d))
    vals.required = (
        if (freq == "hour") hourly.vals.required else
        if (freq == "half hour") 2*hourly.vals.required else
        if (freq == "10 min") 6*hourly.vals.required else stop())
    d = d[, by = .(date, stn),
       if (.N >= vals.required)
          {if (variable.name == "temp.C")
               .(
                    temp.C.mean = mean(value),
                    temp.C.max = max(value),
                    temp.C.min = min(value))
           else if (variable.name == "precipitation.mm")
                .(precipitation.mm.total = sum(value))
           else
                .(value.mean = mean(value))}]
    if ("value.mean" %in% names(d))
        setnames(d, "value.mean", paste0(variable.name, ".mean"))
    d}

combine.dailies = function(l)
    Reduce(x = l, function(x, y)
        merge(x, y, all = T, by = c("stn", "date")))

spanish.month.abbrs = c("ene", "feb", "mar", "abr", "may", "jun", "jul", "ago", "sep", "oct", "nov", "dic")

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

get.ground.raw.simat = function()
  # According to an email exchange Iván Gutiérrez-Avila had with
  # somebody at SIMAT, all times are in UTC-06:00.
  # Info:
  #   Air pressure: http://web.archive.org/web/20190412142523/http://www.aire.cdmx.gob.mx/descargas/datos/excel/Presion-Atmosferica.pdf
  #   Everything else: https://web.archive.org/web/20190412141446/http://www.aire.cdmx.gob.mx/descargas/datos/excel/REDMETxls.pdf
  # There's precipitation data available, but only with weekly samples:
  #   http://148.243.232.112:8080/opendata/redda/concentracion.csv - code "PP"
   {years = year(earliest.date) : (year(Sys.Date()) - 1)
    url.root = "http://148.243.232.112:8080/opendata"

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

    punl(stations, obs)}
get.ground.raw.simat = pairmemo(get.ground.raw.simat, pairmemo.dir)

## ** UNAM

get.ground.raw.unam = function()
   {url.root = "http://www.ruoa.unam.mx/pembu/Estaciones"
    stns = c("CCA", paste0("ENP", 1:9), paste0("CCH", c("A", "N", "O", "S", "V")))
    stations = data.table()
    message("Loading UNAM")
    obs = rbindlist(unlist(rec = F, pblapply(stns, function(stn.nominal)
      # The station is "nominal" in the sense that files for one
      # station can have observations that are stated inside the file
      # to come from a different station.
       {years = (
            (if (stn.nominal == "CCA") 2007 else year(earliest.date)) :
              # CCA observations start on this year.
            (year(Sys.Date()) - 1))
        unlist(rec = F, lapply(years, function(dyear)
           {files = download(
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
                        unname(slurp.zip(path))
                    else
                      # XTRA
                        NULL})
            if (is.null(files))
                return(data.table())

            # Sometimes, files seem to consist of multiple CSV
            # files squished together, so break them apart.
            # Other files are empty, so filter them out.
            files = unlist(lapply(files, function(text)
                if (str_length(text))
                    str_split(text, "\nPrograma de Estaciones")[[1]]
                else
                    character(0)))

            lapply(files, function(text)
               {stopifnot(str_match(text, mlr("^Horario (\\S+)"))[,2]
                    == "UTC-6h")

                stn.name = str_replace_all(
                    str_match(text, mlr("^Estacion (.+?)-UNAM,Ciudad de Mexico"))[,2],
                    "\\s+", "")
                if (stn.name == "CCH0")
                    stn.name = "CCHO"
                latlon = str_match(text, mlr("^Lat ([0-9.]+) N, Lon ([0-9.]+) W"))
                station = data.table(
                    stn.name = stn.name,
                    lon = -as.numeric(latlon[,3]),
                    lat = as.numeric(latlon[,2]))
                stations <<- unique(rbind(stations, station))
                stn = stations[, which(
                    stn.name == station$stn.name &
                    lon == station$lon &
                    lat == station$lat)]

                d = fread(
                    str_sub(text, str_locate(text, mlr("^[0-9]{4}/"))[,1]),
                    col.names = str_split(str_extract(text, mlr("^Fecha_hora,.+")), ",")[[1]],
                    na.strings = "null")

                # Read dates, and drop illegal dates like September 31st.
                d[, date := as.Date(Fecha_hora)]
                d = d[!is.na(date)]

                d = d[, by = date, .(
                    stn = stn,
                    temp.C.mean = if.enough.halfhourly(mean, Temp),
                    temp.C.max = if.enough.halfhourly(max, Temp),
                    temp.C.min = if.enough.halfhourly(min, Temp),
                    wind.speed.sustained.mps.mean = if.enough.halfhourly(mean, Rapidez_v_sostenido),
                    wind.speed.gust.mps.mean = if.enough.halfhourly(mean, Rapidez_rachas),
                    relative.humidity.percent.mean = if.enough.halfhourly(mean, Hum_Rel),
                    pressure.hPa.mean = if.enough.halfhourly(mean, Presion_bar))]
                # Drop rows that are missing on all the data columns.
                vnames = setdiff(colnames(d), c("date", "loc"))
                d[rowSums(!is.na(d[, mget(vnames)])) > 0]})}))})))

    stations[, stn := .I]
    punl(stations, obs)}
get.ground.raw.unam = pairmemo(get.ground.raw.unam, pairmemo.dir)

## ** SMN

## *** Observatories

get.ground.raw.smn.observatories = function()
  # We use the hourly file. This network also provides daily and
  # 15-minute files. However, judging from the ELEMENT-CODE
  # column of these files, they may not have the same variables
  # available.
   {message("Loading SMN observatories")

    na.value = -99999

    # Read in the the giant CSV of hourly observations.
    whole = fread(cmd = paste("unzip -p", shQuote(file.path(
        data.root, "simat-raw", "OBS_2018", "Observatorios_Horarios.zip"))))
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
       {o = file(encoding = "Windows-1252", file.path(data.root,
            "simat-raw", "OBS_2018",
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
    stations = as.data.table(read_excel(file.path(data.root,
        "simat-raw", "OBS_2018", "Claves_Observatorios.xlsx")))
    stations = stations[!is.na(CLAVE), .(
        stn = as.integer(CLAVE),
        lon = -(LONG + LONM/60 + ifelse(is.na(LONS), 0, LONS/(60*60))),
        lat = LATG + LATM/60 + ifelse(is.na(LATS), 0, LATS/(60*60)))]
    
    punl(stations, obs)}
get.ground.raw.smn.observatories = pairmemo(get.ground.raw.smn.observatories, pairmemo.dir)

## *** ESIMES and EMAS

get.ground.raw.smn.esimes = function()
    es.stations(process.es.observations(read.es()))
get.ground.raw.smn.esimes = pairmemo(get.ground.raw.smn.esimes, pairmemo.dir)

get.ground.raw.smn.emas = function()
    es.stations(emas = T,
        process.es.observations(emas = T, n.jobs = 8, read.es(emas = T)))
get.ground.raw.smn.emas = pairmemo(get.ground.raw.smn.emas, pairmemo.dir)

read.es = function(emas = F)
   {fnames = list.files(
       (if (emas)
           file.path(data.root, "simat-emas-csv") else
           file.path(data.root, "simat-raw", "ESIMEs_2018")),
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
    rbindlist(fill = T, Filter(nrow, pblapply(ds, cl = n.jobs, function(d)
       {di <<- di + 1
        #message(di)
        fname = names(ds)[di]

        # Skip empty files.
        if (!nrow(d))
            return(d)

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

        ifcol = function(...)
            for (col in c(...))
                if (col %in% colnames(d))
                   {v = d[[col]]
                    if (is.character(v))
                        v = as.numeric(str_replace(v, fixed(","), "."))
                    return(v)}
        d = d[, .(
            stn,
            date = as.Date(tz = target.tz,
              # Iván Gutiérrez-Avila says all timestamps are in UTC.
              # This is plausible given the observed temperatures.
               {if (is.character(DateTime))
                   as.POSIXct(
                       str_replace_all(DateTime, fixed("."), ""),
                       format = (if (any(str_detect(DateTime, "\\d{8}")))
                           "%Y%m%d %H%M" else
                           paste(
                               (if (any(str_detect(DateTime, "/(?:[2-3][0-9]|1[3-9])/")))
                                   "%m/%d/%Y" else
                                   "%d/%m/%Y"),
                               (if (any(str_detect(DateTime, " P\\.?M")))
                                   "%I:%M:%S %p" else
                                   "%H:%M:%S"))),
                       tz = "UTC")
                else
                   {stopifnot(attr(DateTime, "tz") == "UTC")
                    DateTime}}),
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
            wind.speed.kmph = ifcol("RapViento"),
              # The unit for RapViento is a guess based on the units
              # for other EMAs files.
            wind.speed.sustained.kmph = ifcol("WSK(kph)"), 
            wind.speed.gust.kmph = ifcol("WSMK(kph)"),
            wind.speed.u.mps = ifcol("AvgWSU(m/s)"),
            wind.speed.v.mps = ifcol("AvgWSV(m/s)"),
            pressure.hPa = ifcol("BP(mbar)", "AvgBP(mbar)", "PresBarometric"))]

        combine.dailies(lapply(
            setdiff(colnames(d), c("stn", "date")),
            function(vname)
                daily.summary(
                    d[, .(stn, date, value = get(vname))],
                    (if (round(1 / mean(is.na(d[[vname]]))) == 5)
                      # This variable is only provided on an hourly
                      # basis.
                        "hour" else
                        "10 min"),
                    vname)))})))}

es.stations = function(obs, emas = F)
   {stations = as.data.table(
        if (emas) download(
            "https://smn.cna.gob.mx/tools/DATA/Observando%20el%20Tiempo/EMAS/EMAS.xlsx",
            "emas_stations.xlsx",
            read_excel)
        else download(
            "http://web.archive.org/web/20180809001831id_/http://smn1.conagua.gob.mx/emas/catalogoa.html",
            "esimes_stations.html",
            function(fname)
                readHTMLTable(fname,
                    encoding = "UTF-8", which = 6,
                    header = T, stringsAsFactors = F)))

    coord = function(s)
       {s = str_replace_all(s, fixed(","), ".")
        e = apply(MARGIN = 2, FUN = as.numeric,
            str_match(s, "(\\d+)[°º]\\s*(\\d+)['´]\\s*([0-9.]+)")[, -1])
        e[,1] + e[,2]/60 + e[,3]/(60*60)}
    stations = stations[, .(
        stn = NOMBRE,
        lon = -coord(if (emas) `Longitud W` else `Longitud (Oeste)`),
        lat = coord(if (emas) `Latitud N` else `Latitud (Norte)`))]
    # When stations are duplicated, keep only the first.
    stations = stations[, by = stn, head(.SD, 1)]
    stopifnot(!anyNA(stations))

    l = list(seen = unique(obs$stn), known = stations$stn)
    l$seen.orig = l$seen

    for (x in c("seen", "known"))
        {l[[x]] = str_replace(l[[x]], "\\xd1", "")
         l[[x]] = iconv(tolower(l[[x]]), to = "ASCII//TRANSLIT")

         replacements = c(
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
             "\\bescuela nacional de ciencias biol.?gicas(\\z| ipn)" = "escuela nacional de ciencias biologicas i",
             "\\bencb" = "escuela nacional de ciencias biologicas",
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
             "\\bv carranza" = "venustiano carranza",
             "\\bnochix\\b" = "nochistlan",
             "^presa cuchi" = "presa el cuchillo",
             "(presa )?\\bla cang.*" = "presa cangrejera",
             "\\bjuani\\S*" = "juanico",
             "\\btantak\\b" = "tantaquin",
             "\\btlapa\\b.*" = "tlapa de comonfort")
         for (r in names(replacements))
             l[[x]] = str_replace(l[[x]], r, replacements[[r]])

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
        "bahia de kino ii", "^apan", "^calak(mul)?$", "^matat$",
        "^teziu$")
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

# Wunderground's dates are in America/Mexico_City, so we aren't
# actually using it for now.

get.ground.raw.wunderground = function()
   {db = dbConnect(RSQLite::SQLite(), wunderground.db.path)
    stations = as.data.table(dbGetQuery(db, "select
            stn, lon, lat
        from Stations"))
    obs = as.data.table(dbGetQuery(db, "select
            stn,
            date,
            max_temperature as 'temp.F.max',
            min_temperature as 'temp.F.min',
            temperature as 'temp.F.mean',
            precip_today,
            wind_speed,
            pressure
        from Daily"))
    dbDisconnect(db)
    obs[, date := as.Date(as.character(date), format = "%Y%m%d")]
    stopifnot(!anyNA(obs$date))
    punl(stations, obs)}

## ------------------------------------------------------------
## * Main entry point
## ------------------------------------------------------------

get.ground.raw = function()
   {networks = list(
        simat = get.ground.raw.simat(),
        unam = get.ground.raw.unam(),
        smno = get.ground.raw.smn.observatories(),
        esimes = get.ground.raw.smn.esimes(),
        emas = get.ground.raw.smn.emas())

    # Add network identifiers.
    for (net in names(networks))
       {stations = copy(networks[[net]]$stations)
        stations[, network := net]
        obs = copy(networks[[net]]$obs)
        obs[, stn := paste(stn, net)]
        networks[[net]] = punl(stations, obs)}

    # Combine all networks.
    stations = rbindlist(lapply(networks, function(x)
        x$stations[, .(stn, network, lon, lat)]))
    obs = rbindlist(fill = T, lapply(networks, function(x)
        x$obs))

    # Convert units.
    for (unit in list(
            list(from = "F", to = "C"),
            list(from = "mmHg", to = "hPa"),
            list(from = "kmph", to = "mps")))
       {reg = sprintf("\\.%s\\.(mean|min|max|total)$", unit$from)
        cn = copy(colnames(obs))
        for (col.from in cn)
            if (str_detect(col.from, reg))
               {col.to = str_replace(col.from, reg,
                    sprintf(".%s.\\1", unit$to))
                obs[, (col.to) := ifelse(is.na(get(col.to)),
                    (if (unit$from == "kmph" && unit$to == "mps")
                        conv_multiunit(get(col.from), "km / hr", "m / sec") else
                        conv_unit(get(col.from), unit$from, unit$to)),
                    get(col.to))]
                obs[, (col.from) := NULL]}}

    punl(stations, obs)}
