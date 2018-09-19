# Add all packages.
library(lme4)
library(reshape)
library(foreign) 
library(ggplot2)
library(plyr)
library(data.table)
library(Hmisc)
library(mgcv)
library(gdata)
library(car)
library(dplyr)
library(ggmap)
library(broom)
library(splines)
library(DataCombine)
library(FNN)
library(gstat)

# Functions to join datasets via space and time.
source("./code/r.scripts/geomerge_alpha.r")

# Import clipped grid of Mexico study area.
fullgrid <- fread("./data/work/mexico_grid_ndvi_water_final.csv")
fullgrid$lstid <- paste(fullgrid$long_lst, fullgrid$lat_lst,  sep = "-")
# Original code below. I assume wrong columns were selected for ndviid?
# fullgrid$ndviid<-paste(fullgrid$long_lst_1, fullgrid$lat_lst_1,  sep = "-")
fullgrid$ndviid <- paste(fullgrid$long_ndvi, fullgrid$lat_ndvi, sep = "-")

# Load meteorological station dataset.
Temp <- readRDS("./data/work/all_stations_final.rds")
Temp <- filter(Temp,hi.temp != "NA")
Temp <- as.data.table(Temp)
Temp[, day:= as.Date(strptime(date, "%Y-%m-%d"))]
Temp[, c := as.numeric(format(day, "%Y")) ]

# Load LST data (based on michael dorman R script).
aqua.2013 <- readRDS("./data/RAW/MODIS.AQUA.TERRA.LST.NDVI/stage2/MYD11A1_2013.rds")
a.2013 <- subset(aqua.2013, lstid %in% fullgrid$lstid)
a.2013 <- as.data.table(a.2013)
# Create full LU TS
days <- seq.Date(from = as.Date("2013-01-01"), to = as.Date("2013-12-31"), 1)
# Create date range
days2013 <- data.table(expand.grid(lstid = fullgrid[, unique(lstid)], day = days))
# days2013$lstid <- as.character(days2013$lstid)
# Merge 
setkey(a.2013, lstid, day)
setkey(days2013,lstid, day)
db2013 <- merge(days2013, a.2013, all.x = T)

# Subset grid. Select desired variables/columns.
fullgrid <- select(fullgrid, lstid, elevation, aspectmean, roaddenmean, openplace, ndviid, in_water)

####### Spatial 
# Bring in all spatial components
# Merge
setkey(db2013, lstid)
setkey(fullgrid, lstid)
db2013 <- merge(db2013, fullgrid, all.x = T)  
gc()
head(db2013)

# Ignore this section since % open place is in dataset.
##### import the open places percent , csv in to the database!!!!!
##take note this will be missing in PA areas
# open<-fread("/media/NAS/Uni/Projects/P045_Israel_LST/2.work/open_places.csv")
# open<-select(open,LSTID, PERCENT)
# setnames(open,"LSTID","lstid")
# setnames(open,"PERCENT","open_place_percent")
# setkey(db2013,lstid)
# setkey(open ,lstid)
# db2013 <- merge(db2013,open, all.x = T) 
#############


# Add month.
db2013[, m := as.numeric(format(day, "%m"))]
# Add season.
# 1-winter, 2-spring, 3-summer, 4-autumn
db2013$season <- recode(db2013$m,'1'="1",'2'="1",'3'="2",'4'="2",'5'="2",'6'="3",'7'="3",'8'="3",'9'="4",'10'="4",'11'="4",'12'="1" )
# 1-winter, 2-summer
db2013$seasonSW <- recode(db2013$m,'1'="1",'2'="1",'3'="1",'4'="2",'5'="2",'6'="2",'7'="2",'8'="2",'9'="2",'10'="1",'11'="1",'12'="1")


# Join NDVI to LST.
fin.ndvi <- readRDS("./data/RAW/MODIS.AQUA.TERRA.LST.NDVI/stage2/MYD13A3_2013.rds")
fin.ndvi <- as.data.table(fin.ndvi)
fin.ndvi[, m := as.numeric(format(day, "%m"))]
# fin.ndvi <- filter(fin.ndvi, c==2013)
names(fin.ndvi)[1] <- paste("ndviid")
f.ndvi<-subset(fin.ndvi, ndviid %in% fullgrid$ndviid)

# Add NDVI.
setkey(db2013, ndviid, m)
setkey(f.ndvi, ndviid, m)
db2013 <- merge(db2013, f.ndvi[,list(ndviid,ndvi,m)], all.x = T)
db2013 <- db2013[complete.cases(db2013$lat_lst),]
gc()
summary(db2013)


##### change the date time string  to just  4 digits string)
#Temp$date=as.numeric(Temp$date)
#Temp<-select(Temp,  stn, date,  tempcmax ,tempcmin ,tempcmean,rhmean  , wsmean ,wdmean,   name  ,  itm_e , itm_n, lat, elev)
#Temp$date<- sub("-[[:digit:]]+","",Temp$date)
#Temp$date<- sub("-[[:digit:]]+","",Temp$date)
#Temp$date=as.numeric(Temp$date)


# Relative humidity mean
Temp2013 <- filter(Temp, c==2013)
temp2013tc <- select(Temp2013, stn, r.humidity.mean, lat_stn=Y, long_stn=X, day)
temp2013tc <- na.omit(temp2013tc)
temp2013tc <- as.data.table(temp2013tc)
temp2013tc$stn <- as.character(temp2013tc$stn)

# Spatiotemporal join. 
# Create point matrices for temperature. 
met.m <- makepointsmatrix(temp2013tc, "long_stn", "lat_stn", "stn")
setkey(db2013, lstid)
lu.m <- makepointsmatrix(db2013[db2013[,unique(lstid)], list(long_lst, lat_lst, lstid), mult = "first"], "long_lst", "lat_lst", "lstid")

# Use nearestbyday function to find nearest met station by day with relative humidity mean data.
closestaodse <- nearestbyday(lu.m, met.m, 
                            db2013, temp2013tc[, list(day,r.humidity.mean,stn)], 
                            "lstid", "stn", "meanT", "r.humidity.mean", knearest = 7, maxdistance = 50000)


setkey(db2013, lstid, day)
setkey(closestaodse, lstid, day)
db2013 <- merge(db2013, closestaodse[,list(day, r.humidity.mean, lstid)], all.x = T)

# Wind speed mean
Temp2013 <- filter(Temp,c==2013)
temp2013tc <- select(Temp2013,stn,wind.speed.mean,lat_stn= Y ,long_stn= X,day)
temp2013tc <- na.omit(temp2013tc)
temp2013tc <- as.data.table(temp2013tc)
temp2013tc$stn <- as.character(temp2013tc$stn)

# Spatiotemporal join.
# Matrices for temperature.
met.m <- makepointsmatrix(temp2013tc, "long_stn", "lat_stn", "stn")
setkey(db2013, lstid)
lu.m <- makepointsmatrix(db2013[db2013[,unique(lstid)], list(long_lst, lat_lst, lstid), mult = "first"], "long_lst", "lat_lst", "lstid")

# Use nearestbyday function to find nearest met station by day with wind speed data.
closestaodse <- nearestbyday(lu.m ,met.m , 
                            db2013, temp2013tc[, list(day, wind.speed.mean, stn)], 
                            "lstid", "stn", "meanT", "wind.speed.mean", knearest = 7, maxdistance = 50000)


setkey(db2013, lstid, day)
setkey(closestaodse, lstid, day)
db2013 <- merge(db2013, closestaodse[,list(day, wind.speed.mean, lstid)], all.x = T)

# Bar (barometric pressure?) mean
Temp2013 <- filter(Temp, c==2013)
temp2013tc <- select(Temp2013, stn, bar.mean, lat_stn=Y, long_stn=X, day)
temp2013tc <- na.omit(temp2013tc)
temp2013tc <- as.data.table(temp2013tc)
temp2013tc$stn <- as.character(temp2013tc$stn)

# Spatiotemporal join. 
# Matrices for temperature. 
met.m <- makepointsmatrix(temp2013tc, "long_stn", "lat_stn", "stn")
setkey(db2013, lstid)
lu.m <- makepointsmatrix(db2013[db2013[,unique(lstid)], list(long_lst, lat_lst, lstid), mult = "first"], "long_lst", "lat_lst", "lstid")

# Use nearestbyday function to find nearest met station by day with wind speed data.
closestaodse<- nearestbyday(lu.m ,met.m , 
                            db2013, temp2013tc[, list(day,bar.mean,stn)], 
                            "lstid", "stn", "meanT", "bar.mean", knearest = 7, maxdistance = 50000)


setkey(db2013, lstid, day)
setkey(closestaodse, lstid, day)
db2013 <- merge(db2013, closestaodse[,list(day, bar.mean, lstid)], all.x = T)


# Rain mean.
Temp2013 <- filter(Temp, c==2013)
temp2013tc <- select(Temp2013, stn, rain.mean, lat_stn=Y, long_stn=X, day)
temp2013tc <- na.omit(temp2013tc)
temp2013tc <- as.data.table(temp2013tc)
temp2013tc$stn <- as.character(temp2013tc$stn)

# Spatiotemporal join.
# Matrices for temperature. 
met.m <- makepointsmatrix(temp2013tc, "long_stn", "lat_stn", "stn")
setkey(db2013, lstid)
lu.m <- makepointsmatrix(db2013[db2013[,unique(lstid)], list(long_lst, lat_lst, lstid), mult = "first"], "long_lst", "lat_lst", "lstid")


# Use nearestbyday function to find nearest met station by day with rain data.
closestaodse<- nearestbyday(lu.m ,met.m , 
                            db2013, temp2013tc[, list(day,rain.mean,stn)], 
                            "lstid", "stn", "meanT", "rain.mean", knearest = 7, maxdistance = 50000)


setkey(db2013, lstid, day)
setkey(closestaodse, lstid, day)
db2013 <- merge(db2013, closestaodse[,list(day, rain.mean, lstid)], all.x = T)


Temp <- readRDS("./data/work/all_stations_final.rds")
Temp <- filter(Temp,hi.temp != "NA")
Temp <- as.data.table(Temp)
Temp[, day:=as.Date(strptime(date, "%Y-%m-%d"))]
Temp[, c := as.numeric(format(day, "%Y"))]
Temp2013tc <- as.data.frame(Temp)
Temp2013tc <- select(Temp2013tc, stn, day,hi.temp, low.temp, temp.mean, long_stn=X, lat_stn=Y)



#loop min
for(i in unique(db2013$day)) {
  x<-Temp2013tc[Temp2013tc$day==i, ]
  y= db2013[db2013$day==i, ]
  ##########
  # calculate IDW
  library(gstat)
  #defaults to idw (gstat)
  library(sp)
  coordinates(x) = ~ long_stn + lat_stn
  coordinates(y) = ~ long_lst + lat_lst
  #location statment uneeded since we defined coordinates
  inter = gstat(formula = low.temp ~ 1,  data =x)
  z<-predict(object = inter, newdata = y)
  # head(z)
  db2013$predmin[db2013$day==i] = z$var1.pred
  # spplot(z, "var1.pred", at = 0:330)
}

#loop mean
for(i in unique(db2013$day)) {
  
  x<-Temp2013tc[Temp2013tc$day==i, ]
  y= db2013[db2013$day==i, ]
  ##########
  # calculate IDW
  library(gstat)
  #defaults to idw (gstat)
  library(sp)
  coordinates(x) = ~ long_stn + lat_stn
  coordinates(y) = ~ long_lst + lat_lst
  #location statment uneeded since we defined coordinates
  inter = gstat(formula = temp.mean ~ 1,  data =x)
  z<-predict(object = inter, newdata = y)
  # head(z)
  db2013$predmean[db2013$day==i] = z$var1.pred
  # spplot(z, "var1.pred", at = 0:330)
}


#loop max
for(i in unique(db2013$day)) {
  x<-Temp2013tc[Temp2013tc$day==i, ]
  y= db2013[db2013$day==i, ]
  ##########
  # calculate IDW
  library(gstat)
  #defaults to idw (gstat)
  library(sp)
  coordinates(x) = ~ long_stn + lat_stn
  coordinates(y) = ~ long_lst + lat_lst
  #location statment uneeded since we defined coordinates
  inter = gstat(formula = hi.temp ~ 1,  data =x)
  z<-predict(object = inter, newdata = y)
  # head(z)
  db2013$predmax[db2013$day==i] = z$var1.pred
  # spplot(z, "var1.pred", at = 0:330)
}




# 
# 
# ####-------> mean Ta  for mod 2+3
# 
# #Tempcmean near 7 stations mean
# Temp2013<-filter(Temp,c==2013)
# temp2013tc<-select(Temp2013,stn,tempcmean,lat_stn=  lat ,long_stn=  long,day)
# temp2013tc<-na.omit(temp2013tc)
# temp2013tc$stn<-as.character(temp2013tc$stn)
# 
# 
# #spatio temporal join
# #matrix for temperature 
# met.m <- makepointsmatrix(temp2013tc, "long_stn", "lat_stn", "stn")
# setkey(db2013, lstid)
# lu.m <- makepointsmatrix(db2013[db2013[,unique(lstid)], list(long_lst, lat_lst, lstid), mult = "first"], "long_lst", "lat_lst", "lstid")
# 
# 
# closestmta<- nearestbyday(lu.m ,met.m , 
#                             db2013, temp2013tc[, list(day,tempcmean,stn)], 
#                             "lstid", "stn", "stn.near", "tempcmean", knearest = 7, maxdistance = 50000, nearestmean = T)
# 
# 
# #join to DB
# setkey(db2013,lstid,day)
# setkey(closestmta,lstid,day)
# db2013 <- merge(db2013, closestmta[,list(day,stn.nearmean,lstid)], all.x = T)
# setnames(db2013,"stn.nearmean","meanTa")
# gc()
# 
# #Tempcmaxmean near 7 stations mean
# Temp2013<-filter(Temp,c==2013)
# temp2013tc<-select(Temp2013,stn,tempcmax,lat_stn=  lat ,long_stn=  long,day)
# temp2013tc<-na.omit(temp2013tc)
# temp2013tc$stn<-as.character(temp2013tc$stn)
# 
# 
# #-------> mean Ta  for mod 2+3
# #spatio temporal join
# #matrix for temperature 
# met.m <- makepointsmatrix(temp2013tc, "long_stn", "lat_stn", "stn")
# setkey(db2013, lstid)
# lu.m <- makepointsmatrix(db2013[db2013[,unique(lstid)], list(long_lst, lat_lst, lstid), mult = "first"], "long_lst", "lat_lst", "lstid")
# 
# 
# closestmta<- nearestbyday(lu.m ,met.m , 
#                           db2013, temp2013tc[, list(day,tempcmax,stn)], 
#                           "lstid", "stn", "stn.near", "tempcmax", knearest = 7, maxdistance = 50000, nearestmean = T)
# 
# #join to DB
# setkey(db2013,lstid,day)
# setkey(closestmta,lstid,day)
# db2013 <- merge(db2013, closestmta[,list(day,stn.nearmean,lstid)], all.x = T)
# setnames(db2013,"stn.nearmean","meanTamax")
# gc()
# 
# 
# #Tempcminmean near 7 stations mean
# Temp2013<-filter(Temp,c==2013)
# temp2013tc<-select(Temp2013,stn,tempcmin,lat_stn=  lat ,long_stn=  long,day)
# temp2013tc<-na.omit(temp2013tc)
# temp2013tc$stn<-as.character(temp2013tc$stn)
# 
# 
# #-------> mean Ta  for mod 2+3
# #spatio temporal join
# #matrix for temperature 
# met.m <- makepointsmatrix(temp2013tc, "long_stn", "lat_stn", "stn")
# setkey(db2013, lstid)
# lu.m <- makepointsmatrix(db2013[db2013[,unique(lstid)], list(long_lst, lat_lst, lstid), mult = "first"], "long_lst", "lat_lst", "lstid")
# 
# 
# closestmta<- nearestbyday(lu.m ,met.m , 
#                           db2013, temp2013tc[, list(day,tempcmin,stn)], 
#                           "lstid", "stn", "stn.near", "tempcmin", knearest = 7, maxdistance = 50000, nearestmean = T)
# 
# 
# #join to DB
# setkey(db2013,lstid,day)
# setkey(closestmta,lstid,day)
# db2013 <- merge(db2013, closestmta[,list(day,stn.nearmean,lstid)], all.x = T)
# setnames(db2013,"stn.nearmean","meanTamin")
# gc()
# 
# 

#save
gc()
######## MAKE SURE YOU TOOK OUT STN 
saveRDS(db2013,"./data/outputs/AQUA/2013/c02/MEXICO.mod3.AQ.2013.rds")
gc()


# take out missing night LST >>> mod3 night
# take out missing day LST >>> mod3 night
#create mod 2 file
db2013.m2.day <- db2013[!is.na(d.tempc)]
db2013.m2.night <- db2013[!is.na(n.tempc)]
#rm m3
rm(db2013)
gc()
#save mod2

saveRDS(db2013.m2.day,"./data/outputs/AQUA/2013/c02/MEXICO.mod2.AQ.2013.day.rds")
saveRDS(db2013.m2.night,"./data/outputs/AQUA/2013/c02/MEXICO.mod2.AQ.2013.night.rds")

gc()




########--------->mod1 night
#to fix missing days issues resulting in cartesean error
db2013days <- sort(unique(db2013.m2.night$day))

#Ta import again
Temp <- readRDS("./data/work/all_stations_final.rds")
# Temp <- filter(Temp,hi.temp != "NA")
Temp[, day:=as.Date(strptime(date, "%Y-%m-%d"))]
Temp[, c := as.numeric(format(day, "%Y"))]
Temp2013 <- select(Temp, stn, day, hi.temp, low.temp, temp.mean, c, long_stn=X, lat_stn=Y)
Ta <- filter(Temp2013, c==2013)
Ta <- as.data.table(Ta)
Ta <- select(Ta, stn, day, low.temp, lat_stn, long_stn)

########### join lst to Ta
#create Ta matrix
Ta$stn <- as.character(Ta$stn)
Ta.m <- makepointsmatrix(Ta, "long_stn", "lat_stn", "stn")
#create lst terra matrix
setkey(db2013.m2.night, lstid)
lst.m <- makepointsmatrix(db2013.m2.night[db2013.m2.night[,unique(lstid)], list(long_lst, lat_lst, lstid), mult = "first"], "long_lst", "lat_lst", "lstid")

#run function
closestlst <- nearestbyday(Ta.m, lst.m, 
                           Ta[day %in% db2013days,], db2013.m2.night, 
                           "stn", "lstid", "closest", "n.tempc", knearest = 9, maxdistance = 1500)


#closestlst[,i.stn :=NULL]
closestlst[,closestknn :=NULL]

setkey(Ta, stn, day)
setkey(closestlst, stn, day)
Ta.m1 <- merge(Ta, closestlst, all.x = T)
Ta.m1 <- Ta.m1[!is.na(n.tempc)]
#save mod 1
saveRDS(Ta.m1,"./data/outputs/AQUA/2013/c02/MEXICO.mod1.AQ.2013.night.rds")



########--------->mod1 day
#to fix missing days issues resulting in cartesean error
db2013days <- sort(unique(db2013.m2.day$day))

#Ta import again
Temp <- readRDS("./data/work/all_stations_final.rds")
# Temp<-filter(Temp,hi.temp != "NA")
Temp[, day := as.Date(strptime(date, "%Y-%m-%d"))]
Temp[, c := as.numeric(format(day, "%Y"))]
Temp2013 <- select(Temp, stn, day, hi.temp, low.temp, temp.mean, c,long_stn=X, lat_stn=Y)
Ta <- filter(Temp2013, c==2013)
Ta <- as.data.table(Ta)
Ta <- select(Ta, stn, day, hi.temp, lat_stn,  long_stn)



########### join lst to Ta
#create Ta matrix
Ta$stn <- as.character(Ta$stn)
Ta.m <- makepointsmatrix(Ta, "long_stn", "lat_stn", "stn")
#create lst terra matrix
setkey(db2013.m2.day, lstid)
lst.m <- makepointsmatrix(db2013.m2.day[db2013.m2.day[,unique(lstid)], list(long_lst, lat_lst, lstid), mult = "first"], "long_lst", "lat_lst", "lstid")

#run function
closestlst <- nearestbyday(Ta.m, lst.m, 
                           Ta[day %in% db2013days,], db2013.m2.day, 
                           "stn", "lstid", "closest", "d.tempc", knearest = 9, maxdistance = 1500)


#closestlst[,i.stn :=NULL]
closestlst[,closestknn :=NULL]

setkey(Ta, stn, day)
setkey(closestlst, stn, day)
Ta.m1 <- merge(Ta, closestlst, all.x = T)
Ta.m1 <- Ta.m1[!is.na(d.tempc)]
#save mod 1
saveRDS(Ta.m1,"./data/outputs/AQUA/2013/c02/MEXICO.mod1.AQ.2013.day.rds")



########--------->mod1 mean
#to fix missing days issues resulting in cartesean error
db2013days <- sort(unique(db2013.m2.night$day))

#Ta import again
Temp <- readRDS("./data/work/all_stations_final.rds")
# Temp<-filter(Temp,hi.temp != "NA")
Temp[, day := as.Date(strptime(date, "%Y-%m-%d"))]
Temp[, c := as.numeric(format(day, "%Y"))]
Temp2013 <- select(Temp, stn, day, hi.temp, low.temp, temp.mean, c, long_stn=X, lat_stn=Y)
Ta <- filter(Temp2013, c==2013)
Ta <- as.data.table(Ta)
Ta <- select(Ta, stn, day, temp.mean, lat_stn,  long_stn)

########### join lst to Ta
#create Ta matrix
Ta$stn <- as.character(Ta$stn)
Ta.m <- makepointsmatrix(Ta, "long_stn", "lat_stn", "stn")
#create lst terra matrix
setkey(db2013.m2.night, lstid)
lst.m <- makepointsmatrix(db2013.m2.night[db2013.m2.night[,unique(lstid)], list(long_lst, lat_lst, lstid), mult = "first"], "long_lst", "lat_lst", "lstid")

#run function
closestlst <- nearestbyday(Ta.m, lst.m, 
                           Ta[day %in% db2013days,], db2013.m2.night, 
                           "stn", "lstid", "closest", "n.tempc", knearest = 9, maxdistance = 1500)


#closestlst[,i.stn :=NULL]
closestlst[,closestknn :=NULL]

setkey(Ta, stn, day)
setkey(closestlst, stn, day)
Ta.m1 <- merge(Ta, closestlst, all.x = T)
Ta.m1 <- Ta.m1[!is.na(n.tempc)]
#save mod 1
saveRDS(Ta.m1,"./data/outputs/AQUA/2013/c02/MEXICO.mod1.AQ.2013.mean24.rds")
# ########--------->mod1 mean
# #to fix missing days issues resulting in cartesean error
# db2013days <- sort(unique(db2013.m2.day$day))
# 
# #Ta import again
# Temp<-fread("/media/NAS/Uni/Projects/P045_Israel_LST/2.work/all_stations4.csv")
# # Temp<-filter(Temp,stn != "NA")
# Temp[, day:=as.Date(strptime(date, "%Y-%m-%d"))]
# Temp[, c := as.numeric(format(day, "%Y")) ]
# Ta<-filter(Temp,c==2013)
# Ta<-select(Ta,stn,day,tempcmax, itm_e , itm_n , lat,  long)
# 
# ########### join lst to Ta
# #create Ta matrix
# Ta$stn<-as.character(Ta$stn)
# Ta.m <- makepointsmatrix(Ta, "long", "lat", "stn")
# #create lst terra matrix
# setkey(db2013.m2.day,lstid)
# lst.m <- makepointsmatrix(db2013.m2.day[db2013.m2.day[,unique(lstid)], list(long_lst, lat_lst, lstid), mult = "first"], "long_lst", "lat_lst", "lstid")
# 
# #run function
# closestlst <- nearestbyday(Ta.m, lst.m, 
#                            Ta[day %in% db2013days,], db2013.m2.day, 
#                            "stn", "lstid", "closest", "d.tempc", knearest = 9, maxdistance = 1500)
# 
# 
# #closestlst[,i.stn :=NULL]
# closestlst[,closestknn :=NULL]
# 
# setkey(Ta,stn,day)
# setkey(closestlst,stn,day)
# Ta.m1 <- merge(Ta, closestlst, all.x = T)
# Ta.m1<-Ta.m1[!is.na(d.tempc)]
# #save mod 1
# saveRDS(Ta.m1,"/media/NAS/Uni/Projects/P045_Israel_LST/2.work/mod1.AQ.2013.mean24.rds")



#cleanup
# Warning: you tried to keep "fgrid" which doesn't exist in workspace - nothing was removed
# keep(fgrid, nearestbyday, nearestbydayM1, makepointsmatrix, sure=TRUE) 
keep(fullgrid, nearestbyday, nearestbydayM1, makepointsmatrix, sure=TRUE) 
gc()