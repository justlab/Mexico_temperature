# Load all final temperature results, to combine into one table

library(plyr)
library(data.table)

setwd("J:/PM/Just_Lab/projects/Mexico_temperature")

# Results for max temperature
max2002 <- readRDS("./data/outputs/AQUA/2002/c03/MEXICO.results.2002.AQ.tempmax.rds")
max2003 <- readRDS("./data/outputs/AQUA/2003/c03/MEXICO.results.2003.AQ.tempmax.rds")
max2004 <- readRDS("./data/outputs/AQUA/2004/c03/MEXICO.results.2004.AQ.tempmax.rds")
max2005 <- readRDS("./data/outputs/AQUA/2005/c03/MEXICO.results.2005.AQ.tempmax.rds")
max2006 <- readRDS("./data/outputs/AQUA/2006/c03/MEXICO.results.2006.AQ.tempmax.rds")
max2007 <- readRDS("./data/outputs/AQUA/2007/c03/MEXICO.results.2007.AQ.tempmax.rds")
max2008 <- readRDS("./data/outputs/AQUA/2008/c03/MEXICO.results.2008.AQ.tempmax.rds")
max2009 <- readRDS("./data/outputs/AQUA/2009/c03/MEXICO.results.2009.AQ.tempmax.rds")
max2010 <- readRDS("./data/outputs/AQUA/2010/c03/MEXICO.results.2010.AQ.tempmax.rds")
max2011 <- readRDS("./data/outputs/AQUA/2011/c03/MEXICO.results.2011.AQ.tempmax.rds")
max2012 <- readRDS("./data/outputs/AQUA/2012/c03/MEXICO.results.2012.AQ.tempmax.rds")
max2013 <- readRDS("./data/outputs/AQUA/2013/c03/MEXICO.results.2013.AQ.tempmax.rds")
max2014 <- readRDS("./data/outputs/AQUA/2014/c03/MEXICO.results.2014.AQ.tempmax.rds")
max2015 <- readRDS("./data/outputs/AQUA/2015/c03/MEXICO.results.2015.AQ.tempmax.rds")

# Add column for year
max2002$year <- 2002
max2003$year <- 2003
max2004$year <- 2004
max2005$year <- 2005
max2006$year <- 2006
max2007$year <- 2007
max2008$year <- 2008
max2009$year <- 2009
max2010$year <- 2010
max2011$year <- 2011
max2012$year <- 2012
max2013$year <- 2013
max2014$year <- 2014
max2015$year <- 2015

# Results for mean temperature
mean2002 <- readRDS("./data/outputs/AQUA/2002/c03/MEXICO.results.2002.AQ.tempmean.rds")
mean2003 <- readRDS("./data/outputs/AQUA/2003/c03/MEXICO.results.2003.AQ.tempmean.rds")
mean2004 <- readRDS("./data/outputs/AQUA/2004/c03/MEXICO.results.2004.AQ.tempmean.rds")
mean2005 <- readRDS("./data/outputs/AQUA/2005/c03/MEXICO.results.2005.AQ.tempmean.rds")
mean2006 <- readRDS("./data/outputs/AQUA/2006/c03/MEXICO.results.2006.AQ.tempmean.rds")
mean2007 <- readRDS("./data/outputs/AQUA/2007/c03/MEXICO.results.2007.AQ.tempmean.rds")
mean2008 <- readRDS("./data/outputs/AQUA/2008/c03/MEXICO.results.2008.AQ.tempmean.rds")
mean2009 <- readRDS("./data/outputs/AQUA/2009/c03/MEXICO.results.2009.AQ.tempmean.rds")
mean2010 <- readRDS("./data/outputs/AQUA/2010/c03/MEXICO.results.2010.AQ.tempmean.rds")
mean2011 <- readRDS("./data/outputs/AQUA/2011/c03/MEXICO.results.2011.AQ.tempmean.rds")
mean2012 <- readRDS("./data/outputs/AQUA/2012/c03/MEXICO.results.2012.AQ.tempmean.rds")
mean2013 <- readRDS("./data/outputs/AQUA/2013/c03/MEXICO.results.2013.AQ.tempmean.rds")
mean2014 <- readRDS("./data/outputs/AQUA/2014/c03/MEXICO.results.2014.AQ.tempmean.rds")
mean2015 <- readRDS("./data/outputs/AQUA/2015/c03/MEXICO.results.2015.AQ.tempmean.rds")

# Add column for year
mean2002$year <- 2002
mean2003$year <- 2003
mean2004$year <- 2004
mean2005$year <- 2005
mean2006$year <- 2006
mean2007$year <- 2007
mean2008$year <- 2008
mean2009$year <- 2009
mean2010$year <- 2010
mean2011$year <- 2011
mean2012$year <- 2012
mean2013$year <- 2013
mean2014$year <- 2014
mean2015$year <- 2015

# Results for min temperature
min2002 <- readRDS("./data/outputs/AQUA/2002/c03/MEXICO.results.2002.AQ.tempmin.rds")
min2003 <- readRDS("./data/outputs/AQUA/2003/c03/MEXICO.results.2003.AQ.tempmin.rds")
min2004 <- readRDS("./data/outputs/AQUA/2004/c03/MEXICO.results.2004.AQ.tempmin.rds")
min2005 <- readRDS("./data/outputs/AQUA/2005/c03/MEXICO.results.2005.AQ.tempmin.rds")
min2006 <- readRDS("./data/outputs/AQUA/2006/c03/MEXICO.results.2006.AQ.tempmin.rds")
min2007 <- readRDS("./data/outputs/AQUA/2007/c03/MEXICO.results.2007.AQ.tempmin.rds")
min2008 <- readRDS("./data/outputs/AQUA/2008/c03/MEXICO.results.2008.AQ.tempmin.rds")
min2009 <- readRDS("./data/outputs/AQUA/2009/c03/MEXICO.results.2009.AQ.tempmin.rds")
min2010 <- readRDS("./data/outputs/AQUA/2010/c03/MEXICO.results.2010.AQ.tempmin.rds")
min2011 <- readRDS("./data/outputs/AQUA/2011/c03/MEXICO.results.2011.AQ.tempmin.rds")
min2012 <- readRDS("./data/outputs/AQUA/2012/c03/MEXICO.results.2012.AQ.tempmin.rds")
min2013 <- readRDS("./data/outputs/AQUA/2013/c03/MEXICO.results.2013.AQ.tempmin.rds")
min2014 <- readRDS("./data/outputs/AQUA/2014/c03/MEXICO.results.2014.AQ.tempmin.rds")
min2015 <- readRDS("./data/outputs/AQUA/2015/c03/MEXICO.results.2015.AQ.tempmin.rds")

# Add column for year
min2002$year <- 2002
min2003$year <- 2003
min2004$year <- 2004
min2005$year <- 2005
min2006$year <- 2006
min2007$year <- 2007
min2008$year <- 2008
min2009$year <- 2009
min2010$year <- 2010
min2011$year <- 2011
min2012$year <- 2012
min2013$year <- 2013
min2014$year <- 2014
min2015$year <- 2015

# rbind all results
all_temp_results <- rbind(max2002, max2003, max2004, max2005, max2006, max2007, max2008, 
                          max2009, max2010, max2011, max2012, max2013, max2014, max2015,
                          mean2002, mean2003, mean2004, mean2005, mean2006, mean2007, mean2008, 
                          mean2009, mean2010, mean2011, mean2012, mean2013, mean2014, mean2015,
                          min2002, min2003, min2004, min2005, min2006, min2007, min2008, 
                          min2009, min2010, min2011, min2012, min2013, min2014, min2015)

# Keep only columns with results
all_temp_results <- all_temp_results[, c("year", "type", "m1.R2", "m1.rmspe", "m1.R2.space", "m1.R2.time", "m1.rmspe.space", 
                                         "m1cv.R2", "m1cv.I", "m1cv.Ise", "m1cv.slope", "m1cv.slopese", 
                                         "m1cv.rmspe", "m1cv.R2.space", "m1cv.R2.time", "m1cv.rmspe.space",
                                         "m3.t33", "m3.R2", "m3.rmspe", "m3.R2.space", "m3.R2.time", "m3.rmspe.space",
                                         "m3.I", "m3.Ise", "m3.slope", "m3.slopese")]



# Arrange by year
all_temp_results <- arrange(all_temp_results, desc(year))

# Export to csv
write.csv(all_temp_results, "J:/PM/Just_Lab/projects/Mexico_temperature/figures/r2_results.csv")

# Convert to datatable
all_temp_results_dt <- as.data.table(all_temp_results)

# Calculate mean CV R2 for mod 1
all_temp_results_dt[type=="tempmax", mean(m1cv.R2)]
# [1] 0.8882589
all_temp_results_dt[type=="tempmean", mean(m1cv.R2)]
# [1] 0.9302933
all_temp_results_dt[type=="tempmin", mean(m1cv.R2)]
# [1] 0.8260704

# Calculate mean R2 for mod 3
all_temp_results_dt[type=="tempmax", mean(m3.R2)]
# [1] 0.9399199
all_temp_results_dt[type=="tempmean", mean(m3.R2)]
# [1] 0.9671565
all_temp_results_dt[type=="tempmin", mean(m3.R2)]
# [1] 0.9273144