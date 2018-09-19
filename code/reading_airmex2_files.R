
install.packages("miceadds")
library(miceadds)

setwd("J:/PM/Just_Lab/projects/airmex2/data/monitors/imported")

load.Rdata( filename = "monitoring_data_2004_2016_2016-10-14.RData" , "monitor_data" )
tail(monitor_data) 
