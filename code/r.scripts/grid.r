library(readxl)
library(lme4)
library(reshape)
library(foreign) 
library(ggplot2)
library(plyr)
library(data.table)
library(reshape2)
library(Hmisc)
library(mgcv)
library(gdata)
library(car)
library(dplyr)
library(ggmap)
library(broom)
library(splines)
library(DataCombine)

a<-readRDS("R:\\RAW\\MODIS.AQUA.TERRA.LST.NDVI\\stage1\\MYD11A1_2002.rds")

grid = a[, c("lstid", "lat_lst", "long_lst","day")]

grid [grid$day != "2002-12-31"] = NA
grid<- grid[complete.cases(grid$day),]


write.csv(grid, "R:\\RAW\\grid\\grid.csv")
# 
# x<-unique(grid$lstid,grid$lat_lst,grid$long_lst)
# a <-unique(grid$lstid,incomparables = FALSE)
# b <-unique(grid$lat_lst,incomparables = FALSE)
# c <-unique(grid$long_lst,incomparables = FALSE)
# a<-as.data.table(a)
# b<-as.data.table(b)
# c<-as.data.table(c)
# x <- rbindlist(list(a,b,c ), fill=TRUE)
# colnames(x)[1] = "lstid"
# colnames(x)[2] = "lat_lst"
# colnames(x)[3] = "long_lst"
# summary(x)
