suppressPackageStartupMessages(
   {library(sf)
    library(jsonlite)})

source("../Just_universal/code/pairmemo.R")

data.root = "/data-belle/Mextemp-temporary-testing"
pairmemo.dir = file.path(data.root, "pairmemo")
ground.json.path = file.path(data.root, "ground.json.gz")

crs.lonlat = 4326 # https://epsg.io/4326
crs.mexico.city = 6369 # https://epsg.io/6369
crs.satellite = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"

earliest.date = "2003-01-01"
  # The earliest date we're interested in.
latest.year = 2018L
  # The last year we're interested in.

study.area.buffer.meters = 50e3

pred.area = function()
    st_read(quiet = T, file.path(data.root, "geography", "mxcity_megalopolis"))
      # From https://www.arcgis.com/home/item.html?id=c72bd82a8d6d428bb6914590d6326f7e
pred.area = pairmemo(pred.area, pairmemo.dir, mem = T)

study.area = function()
   {b = st_bbox(st_transform(crs = crs.lonlat,
        st_buffer(pred.area(), study.area.buffer.meters)))
    f = function(z) round(unname(z), 1)
    list(left = f(b$xmin), right = f(b$xmax),
        bottom = f(b$ymin), top = f(b$ymax))}
study.area = pairmemo(study.area, pairmemo.dir, mem = T)

in.study.area = function(lon, lat)
    lon >= study.area()$left & lon <= study.area()$right &
    lat >= study.area()$bottom & lat <= study.area()$top

get.ground = function()
    sapply(simplify = F, fromJSON(gzfile(ground.json.path)),
        as.data.table)
get.ground = pairmemo(get.ground, pairmemo.dir, mem = T)
