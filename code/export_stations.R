# prepare a shapefile of stations to use in QGIS map of figure 1
library(rjson)
library(data.table)
library(sf)

msta = sapply(simplify = F, fromJSON(file = gzfile("/data-belle/Mexico_temperature/ground.json.gz")), as.data.table)
msta = msta$stations
msta[, .N, by = network]
mstaSF = st_as_sf(msta, coords = c("lon", "lat"), crs = 4326)
st_write(mstaSF, "/data-belle/Mexico_temperature/stations_export.shp")
setwd("/data-belle/Mexico_temperature")
zip("/data-belle/Mexico_temperature/stations_shp_export_12-2020.zip", 
    files = list.files(pattern = "stations_export", full.names = TRUE))
# transfer to desktop, then delete exports
unlink(list.files("/data-belle/Mexico_temperature", pattern = "stations_export", full.names = TRUE))
