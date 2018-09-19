#load libraries 
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

# Set working directory.
setwd("/media/Jdrive/PM/Just_Lab/projects/Mexico_temperature")

#function to join on both space and time
source("./code/r.scripts/geomerge_alpha.r")

#import clipped grid
fullgrid<-fread("./data/work/mexico_grid_ndvi_water_final.csv")
fullgrid$lstid<-paste(fullgrid$long_lst, fullgrid$lat_lst,  sep = "-")
fullgrid$ndviid<-paste(fullgrid$long_ndvi, fullgrid$lat_ndvi, sep = "-")
#load met
Temp<-readRDS("./data/work/all_stations_final.rds")
Temp<-filter(Temp,hi.temp != "NA")
Temp<- as.data.table(Temp)
Temp[, day:=as.Date(strptime(date, "%Y-%m-%d"))]
Temp[, c := as.numeric(format(day, "%Y")) ]

#load LST data (BASED on michael dorman R script)
aqua.2002<-readRDS("./data/RAW/MODIS.AQUA.TERRA.LST.NDVI/stage2/MYD11A1_2002.rds")
a.2002<-subset(aqua.2002, lstid %in% fullgrid$lstid)
a.2002<-as.data.table(a.2002)
#create full LU TS
days<-seq.Date(from = as.Date("2002-01-01"), to = as.Date("2002-12-31"), 1)
#create date range
days2002 <-data.table (expand.grid(lstid = fullgrid[, unique(lstid)], day = days))
# days2002$lstid <- as.character(days2002$lstid)
#merge
setkey(a.2002,lstid,day)
setkey(days2002 ,lstid,day)
db2002 <- merge(days2002,a.2002, all.x = T)

#subset fgird, take out unwanted variables/columns
fullgrid<-select(fullgrid,lstid ,elevation , aspectmean  , roaddenmean , openplace,ndviid )

#######spatial 
#bring in all spatial components
    #merge
    setkey(db2002,lstid)
    setkey(fullgrid ,lstid)
    db2002 <- merge(db2002,fullgrid, all.x = T)  
gc()
head(db2002)


##### import the open places percent , csv in to the database!!!!!
##take note this will be missing in PA areas
# open<-fread("/media/NAS/Uni/Projects/P045_Israel_LST/2.work/open_places.csv")
# open<-select(open,LSTID, PERCENT)
# setnames(open,"LSTID","lstid")
# setnames(open,"PERCENT","open_place_percent")
# setkey(db2012,lstid)
# setkey(open ,lstid)
# db2012 <- merge(db2012,open, all.x = T) 
#############



#add month
db2002[, m := as.numeric(format(day, "%m")) ]
#add season
#1-winter, 2-spring,3-summer,4-autum
db2002$season<-recode(db2002$m,"1=1;2=1;3=2;4=2;5=2;6=3;7=3;8=3;9=4;10=4;11=4;12=1" )
#1-winter, 2-summer
db2002$seasonSW<-recode(db2002$m,"1=1;2=1;3=1;4=2;5=2;6=2;7=2;8=2;9=2;10=1;11=1;12=1")


#join NDVI to lst
fin.ndvi<-readRDS("./data/RAW/MODIS.AQUA.TERRA.LST.NDVI/stage2/MYD13A3_2002.rds")
fin.ndvi<-as.data.table(fin.ndvi)
fin.ndvi[, m := as.numeric(format(day, "%m")) ]
# fin.ndvi<-filter(fin.ndvi,c==2002)
names(fin.ndvi)[1]<-paste("ndviid")
f.ndvi<-subset(fin.ndvi, ndviid %in% fullgrid$ndviid)


#add ndvi
setkey(db2002,ndviid,m)
setkey(f.ndvi,ndviid,m)
db2002 <- merge(db2002, f.ndvi[,list(ndviid,ndvi,m)], all.x = T)
db2002<- db2002[complete.cases(db2002$lat_lst),]
gc()
summary(db2002)


##### change the date time string  to just  4 digits string)
#Temp$date=as.numeric(Temp$date)
#Temp<-select(Temp,  stn, date,  tempcmax ,tempcmin ,tempcmean,rhmean  , wsmean ,wdmean,   name  ,  itm_e , itm_n, lat, elev)
#Temp$date<- sub("-[[:digit:]]+","",Temp$date)
#Temp$date<- sub("-[[:digit:]]+","",Temp$date)
#Temp$date=as.numeric(Temp$date)

#rhmean
#rhmean
Temp2002<-filter(Temp,c==2002)
temp2002tc<-select(Temp2002,stn,r.humidity.mean,lat_stn= Y ,long_stn= X,day)
temp2002tc<-na.omit(temp2002tc)
temp2002tc<-as.data.table(temp2002tc)
temp2002tc$stn<-as.character(temp2002tc$stn)

#spatio temporal join
#matrix for temperature 
met.m <- makepointsmatrix(temp2002tc, "long_stn", "lat_stn", "stn")
setkey(db2002, lstid)
lu.m <- makepointsmatrix(db2002[db2002[,unique(lstid)], list(long_lst, lat_lst, lstid), mult = "first"], "long_lst", "lat_lst", "lstid")

#runthescript #rhmean

closestaodse<- nearestbyday(lu.m ,met.m , 
                            db2002, temp2002tc[, list(day,r.humidity.mean,stn)], 
                            "lstid", "stn", "meanT", "r.humidity.mean", knearest = 7, maxdistance = 50000)


setkey(db2002,lstid,day)
setkey(closestaodse,lstid,day)
db2002 <- merge(db2002, closestaodse[,list(day,r.humidity.mean,lstid)], all.x = T)

#wsmean
#wsmean
Temp2002<-filter(Temp,c==2002)
temp2002tc<-select(Temp2002,stn,wind.speed.mean,lat_stn= Y ,long_stn= X,day)
temp2002tc<-na.omit(temp2002tc)
temp2002tc<-as.data.table(temp2002tc)
temp2002tc$stn<-as.character(temp2002tc$stn)

#spatio temporal join
#matrix for temperature 
met.m <- makepointsmatrix(temp2002tc, "long_stn", "lat_stn", "stn")
setkey(db2002, lstid)
lu.m <- makepointsmatrix(db2002[db2002[,unique(lstid)], list(long_lst, lat_lst, lstid), mult = "first"], "long_lst", "lat_lst", "lstid")

#runthescript #wsmean

closestaodse<- nearestbyday(lu.m ,met.m , 
                            db2002, temp2002tc[, list(day,wind.speed.mean,stn)], 
                            "lstid", "stn", "meanT", "wind.speed.mean", knearest = 7, maxdistance = 50000)


setkey(db2002,lstid,day)
setkey(closestaodse,lstid,day)
db2002 <- merge(db2002, closestaodse[,list(day,wind.speed.mean,lstid)], all.x = T)

#bar mean
#bar mean
Temp2002<-filter(Temp,c==2002)
temp2002tc<-select(Temp2002,stn,bar.mean,lat_stn= Y ,long_stn= X,day)
temp2002tc<-na.omit(temp2002tc)
temp2002tc<-as.data.table(temp2002tc)
temp2002tc$stn<-as.character(temp2002tc$stn)

#spatio temporal join
#matrix for temperature 
met.m <- makepointsmatrix(temp2002tc, "long_stn", "lat_stn", "stn")
setkey(db2002, lstid)
lu.m <- makepointsmatrix(db2002[db2002[,unique(lstid)], list(long_lst, lat_lst, lstid), mult = "first"], "long_lst", "lat_lst", "lstid")

#runthescript #bar.mean

closestaodse<- nearestbyday(lu.m ,met.m , 
                            db2002, temp2002tc[, list(day,bar.mean,stn)], 
                            "lstid", "stn", "meanT", "bar.mean", knearest = 7, maxdistance = 50000)


setkey(db2002,lstid,day)
setkey(closestaodse,lstid,day)
db2002 <- merge(db2002, closestaodse[,list(day,bar.mean,lstid)], all.x = T)


#rain mean
#rain mean
Temp2002<-filter(Temp,c==2002)
temp2002tc<-select(Temp2002,stn,rain.mean,lat_stn= Y ,long_stn= X,day)
temp2002tc<-na.omit(temp2002tc)
temp2002tc<-as.data.table(temp2002tc)
temp2002tc$stn<-as.character(temp2002tc$stn)

#spatio temporal join
#matrix for temperature 
met.m <- makepointsmatrix(temp2002tc, "long_stn", "lat_stn", "stn")
setkey(db2002, lstid)
lu.m <- makepointsmatrix(db2002[db2002[,unique(lstid)], list(long_lst, lat_lst, lstid), mult = "first"], "long_lst", "lat_lst", "lstid")



closestaodse<- nearestbyday(lu.m ,met.m , 
                            db2002, temp2002tc[, list(day,rain.mean,stn)], 
                            "lstid", "stn", "meanT", "rain.mean", knearest = 7, maxdistance = 50000)


setkey(db2002,lstid,day)
setkey(closestaodse,lstid,day)
db2002 <- merge(db2002, closestaodse[,list(day,rain.mean,lstid)], all.x = T)


Temp<-readRDS("./data/work/all_stations_final.rds")
Temp<-filter(Temp,hi.temp != "NA")
temp2002tc<-as.data.table(temp2002tc)
Temp[, day:=as.Date(strptime(date, "%Y-%m-%d"))]
Temp[, c := as.numeric(format(day, "%Y")) ]
temp2002tc<-as.data.frame(temp2002tc)
Temp2002<-select(Temp2002,stn,day,hi.temp, low.temp, temp.mean,long_stn= X ,lat_stn= Y)



#loop min
for(i in unique(db2002$day)) {
  x<-Temp2002[Temp2002$day==i, ]
  y= db2002[db2002$day==i, ]
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
  db2002$predmin[db2002$day==i] = z$var1.pred
  # spplot(z, "var1.pred", at = 0:330)
}

#loop mean
for(i in unique(db2002$day)) {
  
  x<-Temp2002[Temp2002$day==i, ]
  y= db2002[db2002$day==i, ]
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
  db2002$predmean[db2002$day==i] = z$var1.pred
  # spplot(z, "var1.pred", at = 0:330)
}


#loop max
for(i in unique(db2002$day)) {
  x<-Temp2002[Temp2002$day==i, ]
  y= db2002[db2002$day==i, ]
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
  db2002$predmax[db2002$day==i] = z$var1.pred
  # spplot(z, "var1.pred", at = 0:330)
}




# 
# 
# ####-------> mean Ta  for mod 2+3
# 
# #Tempcmean near 7 stations mean
# Temp2012<-filter(Temp,c==2012)
# temp2012tc<-select(Temp2012,stn,tempcmean,lat_stn=  lat ,long_stn=  long,day)
# temp2012tc<-na.omit(temp2012tc)
# temp2012tc$stn<-as.character(temp2012tc$stn)
# 
# 
# #spatio temporal join
# #matrix for temperature 
# met.m <- makepointsmatrix(temp2012tc, "long_stn", "lat_stn", "stn")
# setkey(db2012, lstid)
# lu.m <- makepointsmatrix(db2012[db2012[,unique(lstid)], list(long_lst, lat_lst, lstid), mult = "first"], "long_lst", "lat_lst", "lstid")
# 
# 
# closestmta<- nearestbyday(lu.m ,met.m , 
#                             db2012, temp2012tc[, list(day,tempcmean,stn)], 
#                             "lstid", "stn", "stn.near", "tempcmean", knearest = 7, maxdistance = 50000, nearestmean = T)
# 
# 
# #join to DB
# setkey(db2012,lstid,day)
# setkey(closestmta,lstid,day)
# db2012 <- merge(db2012, closestmta[,list(day,stn.nearmean,lstid)], all.x = T)
# setnames(db2012,"stn.nearmean","meanTa")
# gc()
# 
# #Tempcmaxmean near 7 stations mean
# Temp2012<-filter(Temp,c==2012)
# temp2012tc<-select(Temp2012,stn,tempcmax,lat_stn=  lat ,long_stn=  long,day)
# temp2012tc<-na.omit(temp2012tc)
# temp2012tc$stn<-as.character(temp2012tc$stn)
# 
# 
# #-------> mean Ta  for mod 2+3
# #spatio temporal join
# #matrix for temperature 
# met.m <- makepointsmatrix(temp2012tc, "long_stn", "lat_stn", "stn")
# setkey(db2012, lstid)
# lu.m <- makepointsmatrix(db2012[db2012[,unique(lstid)], list(long_lst, lat_lst, lstid), mult = "first"], "long_lst", "lat_lst", "lstid")
# 
# 
# closestmta<- nearestbyday(lu.m ,met.m , 
#                           db2012, temp2012tc[, list(day,tempcmax,stn)], 
#                           "lstid", "stn", "stn.near", "tempcmax", knearest = 7, maxdistance = 50000, nearestmean = T)
# 
# #join to DB
# setkey(db2012,lstid,day)
# setkey(closestmta,lstid,day)
# db2012 <- merge(db2012, closestmta[,list(day,stn.nearmean,lstid)], all.x = T)
# setnames(db2012,"stn.nearmean","meanTamax")
# gc()
# 
# 
# #Tempcminmean near 7 stations mean
# Temp2012<-filter(Temp,c==2012)
# temp2012tc<-select(Temp2012,stn,tempcmin,lat_stn=  lat ,long_stn=  long,day)
# temp2012tc<-na.omit(temp2012tc)
# temp2012tc$stn<-as.character(temp2012tc$stn)
# 
# 
# #-------> mean Ta  for mod 2+3
# #spatio temporal join
# #matrix for temperature 
# met.m <- makepointsmatrix(temp2012tc, "long_stn", "lat_stn", "stn")
# setkey(db2012, lstid)
# lu.m <- makepointsmatrix(db2012[db2012[,unique(lstid)], list(long_lst, lat_lst, lstid), mult = "first"], "long_lst", "lat_lst", "lstid")
# 
# 
# closestmta<- nearestbyday(lu.m ,met.m , 
#                           db2012, temp2012tc[, list(day,tempcmin,stn)], 
#                           "lstid", "stn", "stn.near", "tempcmin", knearest = 7, maxdistance = 50000, nearestmean = T)
# 
# 
# #join to DB
# setkey(db2012,lstid,day)
# setkey(closestmta,lstid,day)
# db2012 <- merge(db2012, closestmta[,list(day,stn.nearmean,lstid)], all.x = T)
# setnames(db2012,"stn.nearmean","meanTamin")
# gc()
# 
# 

#save
gc()
######## MAKE SURE YOU TOOK OUT STN 
saveRDS(db2002,"./data/outputs/AQUA/2002/c02/MEXICO.mod3.AQ.2002.rds")
gc()


# take out missing night LST >>> mod3 night
# take out missing day LST >>> mod3 night
#create mod 2 file
db2002.m2.day <- db2002[!is.na(d.tempc)]
db2002.m2.night <- db2002[!is.na(n.tempc)]
#rm m3
rm(db2002)
gc()
#save mod2

saveRDS(db2002.m2.day,"./data/outputs/AQUA/2002/c02/MEXICO.mod2.AQ.2002.day.rds")
saveRDS(db2002.m2.night,"./data/outputs/AQUA/2002/c02/MEXICO.mod2.AQ.2002.night.rds")

gc()




########--------->mod1 night
#to fix missing days issues resulting in cartesean error
db2002days <- sort(unique(db2002.m2.night$day))

#Ta import again
Temp<-readRDS("./data/work/all_stations_final.rds")
# Temp<-filter(Temp,hi.temp != "NA")
Temp[, day:=as.Date(strptime(date, "%Y-%m-%d"))]
Temp[, c := as.numeric(format(day, "%Y")) ]
Temp2002<-select(Temp,stn,day,hi.temp, low.temp, temp.mean,c,long_stn= X ,lat_stn= Y)
Ta<-filter(Temp2002,c==2002)
Ta<-as.data.table(Ta)
Ta<-select(Ta,stn,day,low.temp, lat_stn,  long_stn)

########### join lst to Ta
#create Ta matrix
Ta$stn<-as.character(Ta$stn)
Ta.m <- makepointsmatrix(Ta, "long_stn", "lat_stn", "stn")
#create lst terra matrix
setkey(db2002.m2.night,lstid)
lst.m <- makepointsmatrix(db2002.m2.night[db2002.m2.night[,unique(lstid)], list(long_lst, lat_lst, lstid), mult = "first"], "long_lst", "lat_lst", "lstid")

#run function
closestlst <- nearestbyday(Ta.m, lst.m, 
                           Ta[day %in% db2002days,], db2002.m2.night, 
                           "stn", "lstid", "closest", "n.tempc", knearest = 9, maxdistance = 1500)


#closestlst[,i.stn :=NULL]
closestlst[,closestknn :=NULL]

setkey(Ta,stn,day)
setkey(closestlst,stn,day)
Ta.m1 <- merge(Ta, closestlst, all.x = T)
Ta.m1<-Ta.m1[!is.na(n.tempc)]
#save mod 1
saveRDS(Ta.m1,"./data/outputs/AQUA/2002/c02/MEXICO.mod1.AQ.2002.night.rds")



########--------->mod1 day
#to fix missing days issues resulting in cartesean error
db2002days <- sort(unique(db2002.m2.day$day))

#Ta import again
Temp<-readRDS("./data/work/all_stations_final.rds")
# Temp<-filter(Temp,hi.temp != "NA")
Temp[, day:=as.Date(strptime(date, "%Y-%m-%d"))]
Temp[, c := as.numeric(format(day, "%Y")) ]
Temp2002<-select(Temp,stn,day,hi.temp, low.temp, temp.mean,c,long_stn= X ,lat_stn= Y)
Ta<-filter(Temp2002,c==2002)
Ta<-as.data.table(Ta)
Ta<-select(Ta,stn,day,hi.temp, lat_stn,  long_stn)



########### join lst to Ta
#create Ta matrix
Ta$stn<-as.character(Ta$stn)
Ta.m <- makepointsmatrix(Ta, "long_stn", "lat_stn", "stn")
#create lst terra matrix
setkey(db2002.m2.day,lstid)
lst.m <- makepointsmatrix(db2002.m2.day[db2002.m2.day[,unique(lstid)], list(long_lst, lat_lst, lstid), mult = "first"], "long_lst", "lat_lst", "lstid")

#run function
closestlst <- nearestbyday(Ta.m, lst.m, 
                           Ta[day %in% db2002days,], db2002.m2.day, 
                           "stn", "lstid", "closest", "d.tempc", knearest = 9, maxdistance = 1500)


#closestlst[,i.stn :=NULL]
closestlst[,closestknn :=NULL]

setkey(Ta,stn,day)
setkey(closestlst,stn,day)
Ta.m1 <- merge(Ta, closestlst, all.x = T)
Ta.m1<-Ta.m1[!is.na(d.tempc)]
#save mod 1
saveRDS(Ta.m1,"./data/outputs/AQUA/2002/c02/MEXICO.mod1.AQ.2002.day.rds")



########--------->mod1 mean
#to fix missing days issues resulting in cartesean error
db2002days <- sort(unique(db2002.m2.night$day))

#Ta import again
Temp<-readRDS("./data/work/all_stations_final.rds")
# Temp<-filter(Temp,hi.temp != "NA")
Temp[, day:=as.Date(strptime(date, "%Y-%m-%d"))]
Temp[, c := as.numeric(format(day, "%Y")) ]
Temp2002<-select(Temp,stn,day,hi.temp, low.temp, temp.mean,c,long_stn= X ,lat_stn= Y)
Ta<-filter(Temp2002,c==2002)
Ta<-as.data.table(Ta)
Ta<-select(Ta,stn,day,temp.mean, lat_stn,  long_stn)

########### join lst to Ta
#create Ta matrix
Ta$stn<-as.character(Ta$stn)
Ta.m <- makepointsmatrix(Ta, "long_stn", "lat_stn", "stn")
#create lst terra matrix
setkey(db2002.m2.night,lstid)
lst.m <- makepointsmatrix(db2002.m2.night[db2002.m2.night[,unique(lstid)], list(long_lst, lat_lst, lstid), mult = "first"], "long_lst", "lat_lst", "lstid")

#run function
closestlst <- nearestbyday(Ta.m, lst.m, 
                           Ta[day %in% db2002days,], db2002.m2.night, 
                           "stn", "lstid", "closest", "n.tempc", knearest = 9, maxdistance = 1500)


#closestlst[,i.stn :=NULL]
closestlst[,closestknn :=NULL]

setkey(Ta,stn,day)
setkey(closestlst,stn,day)
Ta.m1 <- merge(Ta, closestlst, all.x = T)
Ta.m1<-Ta.m1[!is.na(n.tempc)]
#save mod 1
saveRDS(Ta.m1,"./data/outputs/AQUA/2002/c02/MEXICO.mod1.AQ.2002.mean24.rds")
# ########--------->mod1 mean
# #to fix missing days issues resulting in cartesean error
# db2012days <- sort(unique(db2012.m2.day$day))
# 
# #Ta import again
# Temp<-fread("/media/NAS/Uni/Projects/P045_Israel_LST/2.work/all_stations4.csv")
# # Temp<-filter(Temp,stn != "NA")
# Temp[, day:=as.Date(strptime(date, "%Y-%m-%d"))]
# Temp[, c := as.numeric(format(day, "%Y")) ]
# Ta<-filter(Temp,c==2012)
# Ta<-select(Ta,stn,day,tempcmax, itm_e , itm_n , lat,  long)
# 
# ########### join lst to Ta
# #create Ta matrix
# Ta$stn<-as.character(Ta$stn)
# Ta.m <- makepointsmatrix(Ta, "long", "lat", "stn")
# #create lst terra matrix
# setkey(db2012.m2.day,lstid)
# lst.m <- makepointsmatrix(db2012.m2.day[db2012.m2.day[,unique(lstid)], list(long_lst, lat_lst, lstid), mult = "first"], "long_lst", "lat_lst", "lstid")
# 
# #run function
# closestlst <- nearestbyday(Ta.m, lst.m, 
#                            Ta[day %in% db2012days,], db2012.m2.day, 
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
# saveRDS(Ta.m1,"/media/NAS/Uni/Projects/P045_Israel_LST/2.work/mod1.AQ.2012.mean24.rds")



#cleanup
keep(fullgrid,nearestbyday,nearestbydayM1,makepointsmatrix, sure=TRUE) 
gc()
