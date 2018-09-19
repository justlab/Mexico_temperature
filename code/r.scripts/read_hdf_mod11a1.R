library(magrittr)
library(raster)
library(rgdal)
library(gdalUtils)
library(rgeos)
library(reshape2)
library(plyr)

# layers = c(1, 5, 9)

# Set working directory
setwd("/media/michael/Elements/MODIS.AQUA.TERRA.LST.NDVI")
# setwd("~/MODIS.AQUA.TERRA.LST.NDVI")

# START

# Prepare files list
files = list.files(pattern = "^MOD11A1.+\\.hdf$")
files_split = strsplit(files, "\\.")
dates = sapply(files_split, "[", 2)
years = substr(dates, 2, 5)

# Create Extent object for cropping
ext = readOGR("/home/michael/Dropbox/BGU/Adar/processing_mod11A1_mod13A3", "area")
mod_proj = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
ext = spTransform(ext, mod_proj)
# ext = extent(ext) + 15000 # Buffer extent?

###############################################################################################
# STEP 1
###############################################################################################

# # Test
# y = 2001
# d = "A2001185"
# l = 1

for(l in c(1, 5, 9)) {
  
  for(y in 2000:2001) {
    
    result = stack() # Initialize 'result' stack

    for(d in unique(dates[years == y])) {
    
    # Current 4-files group
    current_files = files[dates == d] 
        
      # Empty list for all tiles of given date
      r = list() 
    
      for(k in current_files) {
      
        sds = get_subdatasets(k) # Read current file
        r[[k]] = sds[l] %>% readGDAL %>% raster # Convert to raster
        
      }
      
      tmp = r[[1]]
      
      # If there is more than 1 tile...
      # if(length(current_files) > 1) {
      
      for(k in 2:length(current_files)) {
        
        tmp = mosaic(tmp, r[[k]], fun = "mean") # Mosaic
      
      }
        
      # }
      
      r = crop(tmp, ext) # Crop
      
      names(r) = d # Set layer name
      
      result = stack(result, r) # Add current date mosaic to 'result' stack
      
    }
    
    # Mask
    result = mask(result, mask = ext, progress = "text")
    
    # Read all values into memory
    result = readAll(result)
    
    # Save raster for layer / year
    saveRDS(result, paste0("MOD11A1_", l, "_", y, ".rds"))
    
  }
  
}

###############################################################################################
# STEP 2
###############################################################################################

library(magrittr)
library(raster)
library(rgdal)
library(gdalUtils)
library(rgeos)
library(reshape2)
library(plyr)

# Set working directory- where the RDS files are
setwd("R:\\RAW\\MODIS.AQUA.TERRA.LST.NDVI\\stage1")

for(y in 2000:2015) {
  
  # Gather MODIS images for given year into a list
  d.tempc = readRDS(paste0("MOD11A1_1_", y, ".rds"))
  n.tempc = readRDS(paste0("MOD11A1_5_", y, ".rds"))
  emmisivity = readRDS(paste0("MOD11A1_9_", y, ".rds"))
  
  for(j in 1:nlayers(d.tempc)) {
  
    # List where each item is a 'RasterStack' with ~365 layers
    result = list(d.tempc[[j]], n.tempc[[j]], emmisivity[[j]])
    
    # Convert to 'SpatialPointsDataFrame' and reproject
    result = lapply(result, rasterToPoints, spatial = TRUE)
    result = lapply(result, spTransform, CRSobj = CRS("+proj=longlat +ellps=WGS84 +no_defs"))
    result = lapply(result, as.data.frame)
    
    # Melt
    result = lapply(result, melt, id.vars = c("x", "y"))
    
    # Set variable names
    result[[1]]$var = "d.tempc"
    result[[2]]$var = "n.tempc"
    result[[3]]$var = "emissivity"
    
    # Combine 3 tables
    result = do.call(rbind, result)
    
    # NaN to NA
    result$value[is.nan(result$value)] = NA
    
    # Convert to date
    result$variable = result$variable %>% as.character %>% substr(2, 8) %>% as.Date(format = "%Y%j")
    
    # Cast
    result = dcast(result, ... ~ var, value.var = "value")
    
    # To Celsius
    result$d.tempc = result$d.tempc - 273.15
    result$n.tempc = result$n.tempc - 273.15
    
    # Round coords
    result$x = round(result$x, 3)
    result$y = round(result$y, 3)
    
    # ID
    result$lstid = paste(result$x, result$y, sep = "-")
    
    # Rename
    result = plyr::rename(result, c("variable" = "day", "x" = "long_lst", "y" = "lat_lst"))
    
    # Calculate year
    result$Year = format(result$day, format = "%Y")
    
    # Columns order
    result = result[, c("lstid", "lat_lst", "long_lst", "day", "Year", "d.tempc", "n.tempc", "emissivity")]
    
  }
  
  # write.csv(result, paste0("MOD11A1_", y, ".csv"), row.names = FALSE)
  saveRDS(result, paste0("MOD11A1_", y, ".rds"))
  
}

# Test
# x = grid@data
# names(x) = c("x", "y")
# coordinates(x) = ~ x + y
# x$long_lst = coordinates(x)[, 1]
# x$lat_lst = coordinates(x)[, 2]
# proj4string(x) = "+proj=longlat +ellps=WGS84 +no_defs"
# writeOGR(x, "P:\\2.work", "grid_test2", driver = "ESRI Shapefile")


