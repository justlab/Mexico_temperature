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
                vname.in = "DESCRIPCIO", vname.out = "land.cover"),
            list(url = "https://www.inegi.org.mx/contenido/productos/prod_serv/contenidos/espanol/bvinegi/productos/geografia/tematicas/CLIMAS/702825267568_s.zip",
                fname.out = "climate_types.zip",
                member = "",
                vname.in = "TIPO_C", vname.out = "climate.type",
                crs = 6362)))
       {path = file.path(data.root, "geography", l$fname.out)
        if (!file.exists(path))
           {message("Downloading ", l$fname.out)
            stopifnot(download.file(l$url, path) == 0)}
        dataset = read_sf(paste0("/vsizip/", path, l$member))
        if (is.na(st_crs(dataset)))
            st_crs(dataset) = l$crs
        set(d, j = l$vname.out, value = factor(
            dataset[unlist(st_intersects(
                st_transform(crs = st_crs(dataset), stations.sf),
                dataset)),][[l$vname.in]]))}

    d[network == "wunderground",
        station := sprintf("%03d", as.integer(station))]
    d[order(network, ifelse(network == "wunderground", station, "x"))]}
