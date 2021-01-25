# prepare a shapefile of stations to use in QGIS map of figure 1
source("code/common.R")

library(rjson)
library(data.table)
library(sf)

sta = get.ground()$stations
sta[, .N, by = network]

# to points features
staSF = st_as_sf(sta, coords = c("lon", "lat"), crs = 4326)

# export and ZIP shapefile
out_path = file.path(data.root, "stations_export.shp")
st_write(staSF, out_path)
to_zip = list.files(data.root, pattern = "stations_export", full.names = TRUE)
zip(file.path(data.root, paste0("stations_shp_export_", format(Sys.Date(), "%m-%d-%Y"), ".zip")), 
    files = to_zip)

# delete unzipped shapefile
unlink(to_zip)
