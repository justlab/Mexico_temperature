source("common.R")
source("modeling.R")

station.metadata = function()
   {d = merge(by = "stn",
        stations[, .(
            stn, network,
            station = str_replace(str_replace_all(name, "\\s+", " "),
                " \\S+\\z", ""),
            region = master.grid[stations$mrow, region],
            lon,
            lat,
            elev.m = scales::comma(elev, accuracy = 1))],
        ground[, by = stn, .(
            date.min = min(date), date.max = max(date),
            n.obs = scales::comma(.N))])[, -"stn"]

    # Add land-cover and climate information from separate files.
    stations.sf = st_as_sf(d[, .(lon, lat)],
        crs = crs.lonlat, coords = c(1, 2))
    d[, lon := round(lon, 3)]
    d[, lat := round(lat, 3)]

    for (l in list(

            list(url = "https://www.inegi.org.mx/contenido/productos/prod_serv/contenidos/espanol/bvinegi/productos/geografia/tematicas/uso_suelo/1_250_000/serie_VI/889463598459_s.zip",
                fname.out = "land_cover.zip",
                member = "/conjunto_de_datos/usv250s6_union.shp",
                vname.in = "DESCRIPCIO", vname.out = "land.cover",
                levels = list(
                    "Annual humidity based agriculture" = "agricultura de humedad anual",
                    "Annual irrigation agriculture" = "agricultura de riego anual",
                    "Annual and semi-permanent irrigation agriculture" = "agricultura de riego anual y semipermanente",
                    "Semi-permanent irrigation agriculture" = "agricultura de riego semipermanente",
                    "Annual rainfed agriculture" = "agricultura de temporal anual",
                    "Annual and permanent rainfed agriculture" = "agricultura de temporal anual y permanente",
                    "Human settlements" = "asentamientos humanos",
                    "Oak forest" = "bosque de encino",
                    "Oyamel-fir forest" = "bosque de oyamel",
                    "Pine forest" = "bosque de pino",
                    "Crassicaule shrublands" = "matorral crasicaule",
                    "Sarcocaul shrubland" = "matorral sarcocaule",
                    "Halophilic grassland" = "pastizal halófilo",
                    "Human-induced grassland " = "pastizal inducido",
                    "High mountain meadow" = "pradera de alta montaña",
                    "Secondary (tree type) vegetation of dry broadleaf forest" = "vegetación secundaria arbórea de selva baja caducifolia",
                    "Secondary (bushy type) vegetation of pine-oak forest" = "vegetación secundaria arbustiva de bosque de pino-encino",
                    "Secondary (bushy type) vegetation of dry broadleaf forest" = "vegetación secundaria arbustiva de selva baja caducifolia")),

            list(url = "https://www.inegi.org.mx/contenido/productos/prod_serv/contenidos/espanol/bvinegi/productos/geografia/tematicas/CLIMAS/702825267568_s.zip",
                fname.out = "climate_types.zip",
                member = "",
                vname.in = "TIPO_C", vname.out = "climate.type",
                crs = 6362,
                levels = list(
                    "Warm subhumid" = "cálido subhúmedo",
                    "Cold" = "frío",
                    "Semi-warm subhumid" = "semicálido subhúmedo",
                    "Subhumid semi-cold" = "semifrío subhúmedo",
                    "Semi-dry semi-warm" = "semiseco semicálido",
                    "Temperate semi-dry" = "semiseco templado",
                    "Temperate humid" = "templado húmedo",
                    "Temperate subhumid" = "templado subhúmedo"))))

       {path = file.path(data.root, "geography", l$fname.out)
        if (!file.exists(path))
           {message("Downloading ", l$fname.out)
            stopifnot(download.file(l$url, path) == 0)}
        dataset = read_sf(paste0("/vsizip/", path, l$member))
        if (is.na(st_crs(dataset)))
            st_crs(dataset) = l$crs
        set(d, j = l$vname.out, value = `levels<-`(
            factor(tolower(dataset[unlist(st_intersects(
                st_transform(crs = st_crs(dataset), stations.sf),
                dataset)),][[l$vname.in]])),
            l$levels))
        stopifnot(!anyNA(d[[l$vname.out]]))}

    d[network == "wunderground",
        station := sprintf("%03d", as.integer(station))]
    d[order(network, ifelse(network == "wunderground", station, "x"))]}
