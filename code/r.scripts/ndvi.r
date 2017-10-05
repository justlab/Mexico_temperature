library(magrittr)
library(raster)
library(rgdal)
library(gdalUtils)
library(rgeos)
library(reshape2)
library(plyr)

# Set working directory- where the RDS files are
setwd("R:\\RAW\\MODIS.AQUA.TERRA.LST.NDVI\\stage1")
# setwd("~/Downloads")

for(y in 2002:2002) {
  
  # Gather MODIS images for given year into a list
  ndvi = readRDS(paste0("MYD13A3_1_", y, ".rds"))
  
  final = NULL
  
  for(j in 1:nlayers(ndvi)) {
    
    print(j)
    
    # List where each item is a 'RasterStack' with ~365 layers
    ndvi1 = ndvi[[j]]
    ndvi1 = ndvi1 * 0.0001 * 0.0001
    ndvi1[is.na(ndvi1)] = -9999
    
    # Convert to 'SpatialPointsDataFrame' and reproject
    ndvi1 = rasterToPoints(ndvi1, spatial = TRUE)
    ndvi1 = spTransform(ndvi1, CRSobj = CRS("+proj=longlat +ellps=WGS84 +no_defs"))
    ndvi1 = as.data.frame(ndvi1)
    
    # Melt
    ndvi1 = melt(ndvi1, id.vars = c("x", "y"))
    
    # Convert to date
    ndvi1$variable = ndvi1$variable %>% as.character %>% substr(2, 8) %>% as.Date(format = "%Y%j")
    
    # Back to NA
    ndvi1$value[ndvi1$value == -9999] = NA
    
    # Round coords
    ndvi1$x = round(ndvi1$x, 3)
    ndvi1$y = round(ndvi1$y, 3)
    
    # ID
    ndvi1$lstid = paste(ndvi1$x, ndvi1$y, sep = "-")
    
    # Rename
    ndvi1 = plyr::rename(ndvi1, c("variable" = "day", "x" = "long_lst", "y" = "lat_lst", "value" = "ndvi"))
    
    # Calculate year
    ndvi1$Year = format(ndvi1$day, format = "%Y")
    
    # Columns order
    ndvi1 = ndvi1[, c("lstid", "lat_lst", "long_lst", "day", "Year", "ndvi")]
    
    final = rbind(final, ndvi1)
    
  }
  
  # write.csv(ndvi1, paste0("MOD11A1_", y, ".csv"), row.names = FALSE)
  saveRDS(final, paste0("MYD13A3_", y, ".rds"))
  
}