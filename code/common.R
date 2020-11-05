suppressPackageStartupMessages(
   {library(sf)
    library(jsonlite)
    library(Just.universal)})

data.root = Sys.getenv("JUSTLAB_MEXICO_TEMPERATURE_DATA_ROOT")
stopifnot(file.exists(data.root))

ground.json.path = file.path(data.root, "ground.json.gz")

crs.lonlat = 4326 # https://epsg.io/4326
crs.mexico.city = 6369 # https://epsg.io/6369
crs.satellite = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"

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

earliest.date = "2003-01-01"
  # The earliest date we're interested in.
latest.year = 2019L
  # The last year we're interested in.

study.area.buffer.meters = 50e3

pm = function(...) pairmemo(
    directory = function()
       {dir.create(file.path(data.root, "pairmemo"), showWarnings = F)
        file.path(data.root, "pairmemo")},
    n.frame = 2,
    ...)

pm(mem = T,
pred.area <- function()
    st_read(quiet = T, file.path(data.root, "geography", "mxcity_megalopolis")))
      # From https://www.arcgis.com/home/item.html?id=c72bd82a8d6d428bb6914590d6326f7e

pm(mem = T,
study.area <- function()
   {b = st_bbox(st_transform(crs = crs.lonlat,
        st_buffer(pred.area(), study.area.buffer.meters)))
    f = function(z) round(unname(z), 1)
    list(left = f(b$xmin), right = f(b$xmax),
        bottom = f(b$ymin), top = f(b$ymax))})

in.study.area = function(lon, lat)
    lon >= study.area()$left & lon <= study.area()$right &
    lat >= study.area()$bottom & lat <= study.area()$top

pm(mem = T,
get.ground <- function()
    sapply(simplify = F, fromJSON(gzfile(ground.json.path)),
        as.data.table))
