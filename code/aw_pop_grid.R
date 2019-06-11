# This script calculates area-weighted population from AGEBs that are split across
#   multiple MODIS LST cells. 
# First, it creates a polygon grid from the MODIS point grid. Then it intersects it
#   with the AGEBs and multiplies the AGEB total population by the fraction of area
#   in the smaller intersected feature with the area of the whole AGEB.
# Started by Johnathan Rush on 2019-06-11

library(raster)
library(data.table)
library(sf)
library(sp)
library(mapview)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get the MODIS point grid as a polygon grid ###

# first, get a MODIS HDF to check the raster cell dimensions
# we don't have the MOD11 files for Mexico, we were sent FSTs
# look at a NDVI file instead
ndvi_hdf = list.files("/data-belle/Mexico_temperature/ndvi_c006", pattern = "hdf$", full.names = TRUE)[[1]] # h08v06
ndvi_tile = raster(paste0("HDF4_EOS:EOS_GRID:", ndvi_hdf, ":MOD_Grid_monthly_1km_VI:1 km monthly NDVI"))
xres(ndvi_tile) # 926.6254
yres(ndvi_tile) # 926.6254

# function to create a polygon from centroids, counter-clockwise outer ring
poly_from_centroid <- function(x, y, x_res, y_res){
  mp = matrix(c(x+x_res/2, y+y_res/2, # NE
                x-x_res/2, y+y_res/2, # NW
                x-x_res/2, y-y_res/2, # SW
                x+x_res/2, y-y_res/2, # SE
                x+x_res/2, y+y_res/2 # NE
  ), ncol =2, byrow = T)
  st_polygon(list(mp))
}

# open grid
grid = readRDS(sort(dir(full = T, "/data-belle/Mexico_temperature/pairmemo/get.nonsatellite.data/"))[1])[[1]]
setDT(grid)
grid = grid[x_sinu <= -10.38e6 & x_sinu >= -10.46e6 & y_sinu <= 2.13e6 & y_sinu >= 2.06e6, ]
dim(grid) #   6536  8

# create polygons
grid_polys = mapply(grid[, x_sinu], grid[, y_sinu], xres(ndvi_tile), yres(ndvi_tile), FUN = poly_from_centroid, SIMPLIFY = FALSE)

# sf object with data frame
grid_poly_sf = st_sf(grid[, .(lon, lat, x_sinu, y_sinu, elevation)], 
                     st_sfc(grid_polys, crs = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"))
# rename the very long name of the geometry column
names(grid_poly_sf) <- c(names(grid_poly_sf)[1:(ncol(grid_poly_sf)-1)], "geometry")
st_geometry(grid_poly_sf) <- "geometry"

# Open AGEBs
agebs = read_sf("/home/arferk01/agebmexico2010")
# Reproject to MODIS sinusoidal
agebs_sinu = st_transform(agebs, "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")
# Crop AGEBs to area of interest
bbox <- as(raster::extent(-10.46e6, -10.38e6, 2.06e6, 2.13e6), "SpatialPolygons")
proj4string(bbox) <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
bbox_sf = st_as_sf(bbox)
agebs_bbox = agebs_sinu[st_intersects(agebs_sinu, bbox_sf, sparse = FALSE), ] # keep the agebs that intersect the bounding box polygon
# dim(agebs)      # 56195   198
# dim(agebs_bbox) #   997   198
# reproject cropped AGEBs back to original LCC crs
agebs_bbox = st_transform(agebs_bbox, st_crs(agebs))

# reproject polygon grid to LCC
grid_poly_lcc = st_transform(grid_poly_sf, st_crs(agebs_bbox))
plot(grid_poly_lcc[, "elevation"])

# preview
mapview(agebs_bbox, zcol = "AGEB", legend = FALSE) + mapview(grid_poly_lcc, col.regions = "red", alpha.regions = 0.2, legend = FALSE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Area-weighted population ####

# need a unique ID for each grid cell; use the pairing function on the sinusoidal coordinates in the table
pair<-function(x,y){ 0.5*(x+y)*(x+y+1) +  x }       # pairing function expects positive integers
prepLonSinuM <- function(x){trunc(x+20015110)}      # longitude m, +20015110 makes always positive, as integers
prepLatSinuM <- function(y){trunc(y+10007556)}      # latitude m,  +10007556 makes always positive, as integers
grid_poly_lcc$idLSTpair0 <- pair(grid_poly_lcc$x_sinu, grid_poly_lcc$y_sinu)
grid_poly_lcc

# calculate area of each AGEB
agebs_bbox$area_m2 <- st_area(agebs_bbox)

# intersect AGEBs and polygon grid ~1sec
system.time(ageb_grid <- st_intersection(agebs_bbox, grid_poly_lcc))
dim(agebs_bbox)    #  997  199
dim(grid_poly_lcc) # 6536    7
dim(ageb_grid)     # 3740  205

# get area of intersected AGEBs
ageb_grid$area_intersect <- st_area(ageb_grid)

# area-weighted population of intersected AGEBs
ageb_grid$POBTOT_intersect <- ageb_grid$POBTOT*(ageb_grid$area_intersect/ageb_grid$area_m2)

# sum of population (parts of multiple AGEBs) by grid cell
# can probably do this as a data.frame using dplyr somehow, but I'll do it as a data.table
ageb_gridDT = as.data.table(ageb_grid)
ageb_gridDT[, POBTOT_cell := sum(POBTOT_intersect), by = idLSTpair0]
aw_pop_grid = unique(ageb_gridDT[, .(POBTOT_cell), by = idLSTpair0]) # aw = area-weighted

# copy the sum of the population by grid back to a spatial version of the grid to map it
grid_poly_awpop = merge(grid_poly_lcc, aw_pop_grid, by = "idLSTpair0")

# preview
plot(grid_poly_awpop[, "POBTOT_cell"])
# or
mapview(grid_poly_awpop, zcol = "POBTOT_cell")

