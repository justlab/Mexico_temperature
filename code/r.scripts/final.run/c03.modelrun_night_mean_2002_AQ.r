#add all packages
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

#sourcing
source("./code/r.scripts/CV_splits.r")
source("./code/r.scripts/rmspe.r")

###year 2002
#-------------------->> RES TABLE
res <- matrix(nrow=1, ncol=48)
res <- data.frame(res)
colnames(res) <- c(
  "m1.raw","m1.raw.space","m1.raw.time","m1.time","m1.time.space","m1.time.time","m1.space","m1.space.space","m1.space.time","m1.noaod","m1.noaod.space","m1.noaod.time"
  ,"m1.R2","m1.rmspe","m1.R2.space","m1.R2.time","m1.rmspe.space" #mod1 Full
  ,"m1cv.R2","m1cv.I","m1cv.Ise","m1cv.slope","m1cv.slopese","m1cv.rmspe","m1cv.R2.space","m1cv.R2.time","m1cv.rmspe.space" #mod1 CV
  ,"m1cvloc.R2","m1cvloc.I","m1cvloc.Ise","m1cvloc.slope","m1cvloc.slopese","m1cvloc.rmspe","m1cvloc.R2.space","m1cvloc.R2.time","m1cvloc.rmspe.space"#loc m1
  ,"m2.R2" #mod2
  ,"m3.t31","m3.t33" #mod3 tests
  ,"m3.R2","m3.rmspe","m3.R2.space","m3.R2.time","m3.rmspe.space" #mod3
  ,"m3.I","m3.Ise","m3.slope","m3.slopese")#Extra
res$type <- c("tempmean")

#load data
mod1.n <-readRDS("./data/outputs/AQUA/2002/c02/MEXICO.mod1.AQ.2002.mean24.rds")
summary(mod1.n)

#delete water flags
#kill NA water
mod1.n<-filter(mod1.n,ndvi > 0)
mod1.n<-filter(mod1.n,in_water == 0)
mod1.n<-filter(mod1.n,bar.mean != "NA")
mod1.n<- as.data.table(mod1.n)
# mod1.n<-filter(mod1.n,!is.na(open_place_percent))

#####min temperature

###base linear model night
# l1.raw.formula <- as.formula(tempcmax ~ n.tempc+ndvi+ELEVATION+ASPECT+DENS_POP+roadden+WATER_DIST+rhmean+wsmean+open_place_percent )
# out.l1<-lm(l1.raw.formula,data=mod1.n)
# summary(out.l1)
# mod1$pred.m1 <- predict(out.l1)
# print(summary(lm(tempcmax~pred.m1,data=mod1))$r.squared)


m1.formula <- as.formula(temp.mean ~ n.tempc+ndvi+elevation+aspectmean+roaddenmean+r.humidity.mean+bar.mean+rain.mean+wind.speed.mean+openplace +(1+n.tempc|day))
####
m1_sc <- lmer(m1.formula,data=mod1.n)
summary(m1_sc)
mod1.n$pred.m1 <- predict(m1_sc)
res[res$type=="tempmean", 'm1.R2'] <- print(summary(lm(temp.mean~pred.m1,data=mod1.n))$r.squared)
#RMSPE
res[res$type=="tempmean", 'm1.rmspe'] <- print(rmse(residuals(m1_sc)))

#spatial
spatialall<-mod1.n %>%
  group_by(stn) %>%
  dplyr::summarise(barpm = mean(temp.mean, na.rm=TRUE), barpred = mean(pred.m1, na.rm=TRUE)) 
m1.fit.all.s <- lm(barpm ~ barpred, data=spatialall)
res[res$type=="tempmean", 'm1.R2.space'] <-print(summary(lm(barpm ~ barpred, data=spatialall))$r.squared)
res[res$type=="tempmean", 'm1.rmspe.space'] <- print(rmse(residuals(m1.fit.all.s)))


#temporal
#temporal (take out daily PM from yearly mean)
tempoall<-left_join(mod1.n,spatialall)
tempoall$delpm <-tempoall$temp.mean-tempoall$barpm
tempoall$delpred <-tempoall$pred.m1-tempoall$barpred
mod_temporal <- lm(delpm ~ delpred, data=tempoall)
res[res$type=="tempmean", 'm1.R2.time']<- print(summary(lm(delpm ~ delpred, data=tempoall))$r.squared)


#save
saveRDS(mod1.n,"./data/outputs/AQUA/2002/c03/MEXICO.mod1.2002.AQ.night.mean.predm1.rds")
#save results
saveRDS(res,"./data/outputs/AQUA/2002/c03/MEXICO.results.2002.AQ.tempmean.rds")




#---------------->>>> CV
#s1
splits_s1 <- splitdf(mod1.n)
test_s1 <- splits_s1$testset
train_s1 <- splits_s1$trainset
out_train_s1 <- lmer(m1.formula,data =  train_s1 )
test_s1$pred.m1.cv <- predict(object=out_train_s1 ,newdata=test_s1,allow.new.levels=TRUE,re.form=NULL )
test_s1$iter<-"s1"
#s2
splits_s2 <- splitdf(mod1.n)
test_s2 <- splits_s2$testset
train_s2 <- splits_s2$trainset
out_train_s2 <- lmer(m1.formula,data =  train_s2 )
test_s2$pred.m1.cv <- predict(object=out_train_s2 ,newdata=test_s2,allow.new.levels=TRUE,re.form=NULL )
test_s2$iter<-"s2"
#s3
splits_s3 <- splitdf(mod1.n)
test_s3 <- splits_s3$testset
train_s3 <- splits_s3$trainset
out_train_s3 <- lmer(m1.formula,data =  train_s3 )
test_s3$pred.m1.cv <- predict(object=out_train_s3 ,newdata=test_s3,allow.new.levels=TRUE,re.form=NULL )
test_s3$iter<-"s3"
#s4
splits_s4 <- splitdf(mod1.n)
test_s4 <- splits_s4$testset
train_s4 <- splits_s4$trainset
out_train_s4 <- lmer(m1.formula,data =  train_s4 )
test_s4$pred.m1.cv <- predict(object=out_train_s4 ,newdata=test_s4,allow.new.levels=TRUE,re.form=NULL )
test_s4$iter<-"s4"
#s5
splits_s5 <- splitdf(mod1.n)
test_s5 <- splits_s5$testset
train_s5 <- splits_s5$trainset
out_train_s5 <- lmer(m1.formula,data =  train_s5 )
test_s5$pred.m1.cv <- predict(object=out_train_s5 ,newdata=test_s5,allow.new.levels=TRUE,re.form=NULL )
test_s5$iter<-"s5"
#s6
splits_s6 <- splitdf(mod1.n)
test_s6 <- splits_s6$testset
train_s6 <- splits_s6$trainset
out_train_s6 <- lmer(m1.formula,data =  train_s6 )
test_s6$pred.m1.cv <- predict(object=out_train_s6 ,newdata=test_s6,allow.new.levels=TRUE,re.form=NULL )
test_s6$iter<-"s6"
#s7
splits_s7 <- splitdf(mod1.n)
test_s7 <- splits_s7$testset
train_s7 <- splits_s7$trainset
out_train_s7 <- lmer(m1.formula,data =  train_s7 )
test_s7$pred.m1.cv <- predict(object=out_train_s7 ,newdata=test_s7,allow.new.levels=TRUE,re.form=NULL )
test_s7$iter<-"s7"
#s8
splits_s8 <- splitdf(mod1.n)
test_s8 <- splits_s8$testset
train_s8 <- splits_s8$trainset
out_train_s8 <- lmer(m1.formula,data =  train_s8 )
test_s8$pred.m1.cv <- predict(object=out_train_s8 ,newdata=test_s8,allow.new.levels=TRUE,re.form=NULL )
test_s8$iter<-"s8"
#s9
splits_s9 <- splitdf(mod1.n)
test_s9 <- splits_s9$testset
train_s9 <- splits_s9$trainset
out_train_s9 <- lmer(m1.formula,data =  train_s9 )
test_s9$pred.m1.cv <- predict(object=out_train_s9 ,newdata=test_s9,allow.new.levels=TRUE,re.form=NULL )
test_s9$iter<-"s9"
#s10
splits_s10 <- splitdf(mod1.n)
test_s10 <- splits_s10$testset
train_s10 <- splits_s10$trainset
out_train_s10 <- lmer(m1.formula,data =  train_s10 )
test_s10$pred.m1.cv <- predict(object=out_train_s10 ,newdata=test_s10,allow.new.levels=TRUE,re.form=NULL )
test_s10$iter<-"s10"

#BIND 1 dataset
mod1.n.cv<- data.table(rbind(test_s1,test_s2,test_s3,test_s4,test_s5,test_s6,test_s7,test_s8,test_s9, test_s10))
#save
#saveRDS(mod1.n.cv,"./data/outputs/AQUA/2002/c03/mod1.n.AQ.2002.tempmean.CV.rds")
# cleanup (remove from WS) objects from CV
rm(list = ls(pattern = "train_|test_"))
#table updates
m1.fit.all.cv<-lm(temp.mean~pred.m1.cv,data=mod1.n.cv)
res[res$type=="tempmean", 'm1cv.R2'] <- print(summary(lm(temp.mean~pred.m1.cv,data=mod1.n.cv))$r.squared)
res[res$type=="tempmean", 'm1cv.I'] <-print(summary(lm(temp.mean~pred.m1.cv,data=mod1.n.cv))$coef[1,1])
res[res$type=="tempmean", 'm1cv.Ise'] <-print(summary(lm(temp.mean~pred.m1.cv,data=mod1.n.cv))$coef[1,2])
res[res$type=="tempmean", 'm1cv.slope'] <-print(summary(lm(temp.mean~pred.m1.cv,data=mod1.n.cv))$coef[2,1])
res[res$type=="tempmean", 'm1cv.slopese'] <-print(summary(lm(temp.mean~pred.m1.cv,data=mod1.n.cv))$coef[2,2])
#RMSPE
res[res$type=="tempmean", 'm1cv.rmspe'] <- print(rmse(residuals(m1.fit.all.cv)))


#spatial
spatialall.cv<-mod1.n.cv %>%
  group_by(stn) %>%
  summarise(barpm = mean(temp.mean, na.rm=TRUE), barpred = mean(pred.m1, na.rm=TRUE)) 
m1.fit.all.cv.s <- lm(barpm ~ barpred, data=spatialall.cv)
res[res$type=="tempmean", 'm1cv.R2.space'] <-  print(summary(lm(barpm ~ barpred, data=spatialall.cv))$r.squared)
res[res$type=="tempmean", 'm1cv.rmspe.space'] <- print(rmse(residuals(m1.fit.all.cv.s)))

#temporal
tempoall.cv<-left_join(mod1.n.cv,spatialall.cv)
tempoall.cv$delpm <-tempoall.cv$temp.mean-tempoall.cv$barpm
tempoall.cv$delpred <-tempoall.cv$pred.m1.cv-tempoall.cv$barpred
mod_temporal.cv <- lm(delpm ~ delpred, data=tempoall.cv)
res[res$type=="tempmean", 'm1cv.R2.time'] <-  print(summary(lm(delpm ~ delpred, data=tempoall.cv))$r.squared)


#save
saveRDS(mod1.n.cv,"./data/outputs/AQUA/2002/c03/MEXICO.mod1.2002.AQ.night.mean.predm1.CV.rds")
#save res   
saveRDS(res,"./data/outputs/AQUA/2002/c03/MEXICO.results.2002.AQ.tempmean.rds")


### mod 2 (around 2-4 h)

mod2.n <- readRDS("./data/outputs/AQUA/2002/c02/MEXICO.mod2.AQ.2002.night.rds")
summary(mod2.n)

#delete water flags
mod2.n<-filter(mod2.n,ndvi > 0)
mod2.n<-filter(mod2.n,in_water == 0)
#kill NA water
mod2.n<-filter(mod2.n,bar.mean != "NA")
mod2.n<- as.data.table(mod2.n)
# mod2.n<-filter(mod2.n,!is.na(open_place_percent))
# mod2.n<-filter(mod2.n,!is.na(rhmean))
# mod2.n<-filter(mod2.n,!is.na(wsmean))

mod2.n$pred.m2<-predict(object=m1_sc,newdata=mod2.n,allow.new.levels=TRUE,re.form=NULL)
mod2.n<-filter(mod2.n,!is.na(pred.m2))
mod2.n<- as.data.table(mod2.n)
gc()
setkey(mod2.n,day, lstid)
mod2.n<-mod2.n[!is.na(predmean)]
mod2.n$m <- as.numeric(format(mod2.n$day, "%m")) 
mod2.n[, bimon := (m + 1) %/% 2]
summary(mod2.n$pred.m2)
gc()
mod2.n <- select(mod2.n,day,lstid,m,predmean,long_lst,lat_lst,bimon,pred.m2,n.tempc) ########################################################
saveRDS(mod2.n,"./data/outputs/AQUA/2002/c03/MEXICO.mod2.2002.AQ.night.mean.predm2.rds")
keep(mod2.n,res,rmse,splitdf, sure=TRUE) 
gc()


#aggregate data
#check spatial patterns by plotting a map in mod2
out <-mod2.n %>%
  group_by(lstid) %>%
  summarise(x=mean(long_lst, na.rm=TRUE), y =mean(lat_lst, na.rm=TRUE), pred.m2=mean(pred.m2, na.rm=TRUE)  )
out<-na.omit(out)
write.csv(out,"./data/outputs/AQUA/2002/c03/MEXICO.mod2.2002.AQ.mean.map.csv")

#library(ggmap)

#map <- get_map(location = 'israel', maptype = "roadmap", zoom = 7)
#map <- get_map(location = 'israel', source="stamen", maptype = "watercolor" ,zoom = 7)
#str(MxC_Map_df) #get info on data
#P4 <- ggmap(map, darken = c(0.5, "white"))
#P4 + geom_point(data = out ,aes(x, y,  color = pred.m2))+ ggtitle("predmod2")+scale_colour_gradient2( low = "green", mid = "yellow", high = "red",midpoint = 15, space = "rgb", na.value = "white", guide = "colourbar")

# 
# #mod 3 (5-8 h)
# mod3 <- readRDS("/media/NAS/Uni/Projects/P045_Israel_LST/2.work/mod2.AQ.2004.night.rds")
# mod3x<-select(mod3,meanTamax ,long_lst, lat_lst,night,lstid)
# mod3x<-filter(mod3x,!is.na(meanTamax) )
# 
# mod3x$m <- as.numeric(format(mod3x$night, "%m")) 
# mod3x[, bimon := (m + 1) %/% 2]
# 
# ##test bam
# x <- bam(pred.m2 ~ meanTamax + te(long_lst, lat_lst, by = bimon), data = mod2.n)
# mod3x$pred.m3<- predict.bam(x,mod3x)
# 
# 
# mod1 <-readRDS("/media/NAS/Uni/Projects/P045_Israel_LST/2.work/mod1.min.AQ.2004.tempmean.predm1.rds")
# mod1$lstid<-paste(mod1$long_lst,mod1$lat_lst,sep="-")
# mod1<-mod1[,c("lstid","night","tempcmax","pred.m1","stn"),with=FALSE]
# #R2.m3
# setkey(mod3x,night,lstid)
# setkey(mod1,night,lstid)
# mod1 <- merge(mod1,mod3x[, list(night,lstid,pred.m3)], all.x = T)
# m3.fit.all<- summary(lm(tempcmax~pred.m3,data=mod1))
# print(summary(lm(tempcmax~pred.m3,data=mod1))$r.squared)   



#run lme regression, this *should* include the thin plate spline yet will not run (computational limitations) thus we break it down into 2 components  
summary(mod2.n)
m2.smooth = lme(pred.m2 ~ predmean,random = list(lstid= ~1 + predmean),control=lmeControl(opt = "optim"), data= mod2.n )
#correlate to see everything from mod2.n and the mpm works
mod2.n$pred.t31<-predict(m2.smooth)
mod2.n$resid<-residuals(m2.smooth)
#check R2 
print(summary(lm(pred.m2~pred.t31,data=mod2.n))$r.squared)


#split the files to the separate bi monthly data sets (using dplyr syntax)
# Tall_bimon1 <- filter(mod2.n ,bimon == "1")
# Tall_bimon2 <- filter(mod2.n ,bimon == "2")
# Tall_bimon3 <- filter(mod2.n ,bimon == "3")
Tall_bimon4 <- filter(mod2.n ,bimon == "4")
Tall_bimon5 <- filter(mod2.n ,bimon == "5")
Tall_bimon6 <- filter(mod2.n ,bimon == "6")

#run the separate splines (smooth) for x and y for each bimon
# fit2_1 <- gam(resid ~ s(long_lst,lat_lst),  data= Tall_bimon1 )
# fit2_2 <- gam(resid ~ s(long_lst,lat_lst),  data= Tall_bimon2 )
# fit2_3 <- gam(resid ~ s(long_lst,lat_lst),  data= Tall_bimon3 )
fit2_4 <- gam(resid ~ s(long_lst,lat_lst),  data= Tall_bimon4 )
fit2_5 <- gam(resid ~ s(long_lst,lat_lst),  data= Tall_bimon5 )
fit2_6 <- gam(resid ~ s(long_lst,lat_lst),  data= Tall_bimon6 )

#get the predicted-fitted 
# Xpred_1 <- (Tall_bimon1$pred.t31 - fit2_1$fitted)
# Xpred_2 <- (Tall_bimon2$pred.t31 - fit2_2$fitted)
# Xpred_3 <- (Tall_bimon3$pred.t31 - fit2_3$fitted)
Xpred_4 <- (Tall_bimon4$pred.t31 - fit2_4$fitted)
Xpred_5 <- (Tall_bimon5$pred.t31 - fit2_5$fitted)
Xpred_6 <- (Tall_bimon6$pred.t31 - fit2_6$fitted)

#remerge to 1 file
mod2.n$pred.m2.int <- c(Xpred_4, Xpred_5, Xpred_6)
#this is important so that its sorted as in the first gamm
setkey(mod2.n,day, lstid)

#rerun the lme on the predictions including the spatial spline (smooth)
Final_pred_all <- lme(pred.m2.int ~ predmean ,random = list(lstid= ~1 + predmean ),control=lmeControl(opt = "optim"),data= mod2.n  )
mod2.n$pred.t33 <-predict(Final_pred_all)
#check correlations
res[res$type=="tempmean", 'm3.t33'] <- print(summary(lm(pred.m2 ~ pred.t33,data=mod2.n))$r.squared) 



#mod 3 (5-8 h)
mod3 <- readRDS("./data/outputs/AQUA/2002/c02/MEXICO.mod3.AQ.2002.rds")
#delete water flags
mod3<-filter(mod3,ndvi > 0)
mod3<-filter(mod3,in_water == 0)

#kill NA
mod3<-filter(mod3,!is.na(r.humidity.mean))
mod3<-filter(mod3,!is.na(wind.speed.mean))
mod3<-filter(mod3,!is.na(bar.mean))
mod3<-filter(mod3,!is.na(rain.mean))
mod3<- as.data.table(mod3)

mod3[, m := as.numeric(format(day, "%m")) ]
mod3 <- select(mod3,day,lstid,m,predmean,long_lst,lat_lst)
mod3[, bimon := (m + 1) %/% 2]
setkey(mod3,day, lstid)
mod3<-mod3[!is.na(predmean)]
summary(mod3)
#generate m.3 mix model  predictions 
mod3$pred.m3.mix <-  predict(Final_pred_all,mod3)

#create unique grid
ugrid <-mod3 %>%
  group_by(lstid) %>%
  summarise(long_lst = mean(long_lst, na.rm=TRUE),  lat_lst = mean(lat_lst, na.rm=TRUE)) 

#### PREDICT Gam part
#split back into bimons to include the gam prediction in final prediction        
mod3_bimon1 <- mod3[bimon == 1, ]
mod3_bimon2 <- mod3[bimon == 2, ]
mod3_bimon3 <- mod3[bimon == 3, ]
mod3_bimon4 <- mod3[bimon == 4, ]
mod3_bimon5 <- mod3[bimon == 5, ]
mod3_bimon6 <- mod3[bimon == 6, ]


#addin unique grid to each bimon           
uniq_gid_bimon1 <- ugrid
uniq_gid_bimon2 <- ugrid
uniq_gid_bimon3 <- ugrid
uniq_gid_bimon4 <- ugrid
uniq_gid_bimon5 <- ugrid
uniq_gid_bimon6 <- ugrid

#get predictions for Bimon residuals
# uniq_gid_bimon1$gpred <- predict.gam(fit2_1,uniq_gid_bimon1)
# uniq_gid_bimon2$gpred <- predict.gam(fit2_2,uniq_gid_bimon2)
# uniq_gid_bimon3$gpred <- predict.gam(fit2_3,uniq_gid_bimon3)
uniq_gid_bimon4$gpred <- predict.gam(fit2_4,uniq_gid_bimon4)
uniq_gid_bimon5$gpred <- predict.gam(fit2_5,uniq_gid_bimon5)
uniq_gid_bimon6$gpred <- predict.gam(fit2_6,uniq_gid_bimon6)

#change bimon to data.table
#uniq_gid_bimon1<- as.data.table(uniq_gid_bimon1)
#uniq_gid_bimon2<- as.data.table(uniq_gid_bimon2)
#uniq_gid_bimon3<- as.data.table(uniq_gid_bimon3)
uniq_gid_bimon4<- as.data.table(uniq_gid_bimon4)
uniq_gid_bimon5<- as.data.table(uniq_gid_bimon5)
uniq_gid_bimon6<- as.data.table(uniq_gid_bimon6)

#merge things back togheter
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> merges
# setkey(uniq_gid_bimon1,lstid)
# setkey(mod3_bimon1,lstid)
# mod3_bimon1 <- merge(mod3_bimon1, uniq_gid_bimon1[,list(lstid,gpred)], all.x = T)
# setkey(uniq_gid_bimon2,lstid)
# setkey(mod3_bimon2,lstid)
# mod3_bimon2 <- merge(mod3_bimon2, uniq_gid_bimon2[,list(lstid,gpred)], all.x = T)
# setkey(uniq_gid_bimon3,lstid)
# setkey(mod3_bimon3,lstid)
# mod3_bimon3 <- merge(mod3_bimon3, uniq_gid_bimon3[,list(lstid,gpred)], all.x = T)
setkey(uniq_gid_bimon4,lstid)
setkey(mod3_bimon4,lstid)
mod3_bimon4 <- merge(mod3_bimon4, uniq_gid_bimon4[,list(lstid,gpred)], all.x = T)
setkey(uniq_gid_bimon5,lstid)
setkey(mod3_bimon5,lstid)
mod3_bimon5 <- merge(mod3_bimon5, uniq_gid_bimon5[,list(lstid,gpred)], all.x = T)
setkey(uniq_gid_bimon6,lstid)
setkey(mod3_bimon6,lstid)
mod3_bimon6 <- merge(mod3_bimon6, uniq_gid_bimon6[,list(lstid,gpred)], all.x = T)

#reattach all parts        
mod3 <- rbind(mod3_bimon4,mod3_bimon5,mod3_bimon6)
# create pred.m3    
mod3$pred.m3 <-mod3$pred.m3.mix+mod3$gpred
# hist(mod3$pred.m3)
summary(mod3$pred.m3)
saveRDS(mod3,"./data/outputs/AQUA/2002/c03/MEXICO.mod3.2002.AQ.night.mean.predm3.rds")
keep(mod3,res,rmse, sure=TRUE) 
gc()



#calculate stage 3 R2- CV ten folds approach will take 6 weeks...we don't currently do CV for stage 3.

mod1 <-readRDS("./data/outputs/AQUA/2002/c03/MEXICO.mod1.2002.AQ.night.mean.predm1.rds")
mod1$lstid<-paste(mod1$long_lst,mod1$lat_lst,sep="-")
mod1<-mod1[,c("lstid","day","temp.mean","stn","pred.m1"),with=FALSE]
#R2.m3
setkey(mod3,day,lstid)
setkey(mod1,day,lstid)
mod1 <- merge(mod1,mod3[, list(day,lstid,pred.m3)], all.x = T)
m3.fit.all<- summary(lm(temp.mean~pred.m3,data=mod1))
res[res$type=="tempmean", 'm3.R2'] <- print(summary(lm(temp.mean~pred.m3,data=mod1))$r.squared)    
res[res$type=="tempmean", 'm3.I'] <-print(summary(lm(temp.mean~pred.m3,data=mod1))$coef[1,1])
res[res$type=="tempmean", 'm3.Ise'] <-print(summary(lm(temp.mean~pred.m3,data=mod1))$coef[1,2])
res[res$type=="tempmean", 'm3.slope'] <-print(summary(lm(temp.mean~pred.m3,data=mod1))$coef[2,1])
res[res$type=="tempmean", 'm3.slopese'] <-print(summary(lm(temp.mean~pred.m3,data=mod1))$coef[2,2])
#RMSPE
res[res$type=="tempmean", 'm3.rmspe'] <- print(rmse(residuals(m3.fit.all)))


#spatial
###to check
spatialall<-mod1 %>%
  group_by(stn) %>%
  summarise(barpm = mean(temp.mean, na.rm=TRUE), barpred = mean(pred.m3, na.rm=TRUE)) 
m1.fit.all.spat<- lm(barpm ~ barpred, data=spatialall)
res[res$type=="tempmean", 'm3.R2.space'] <-  print(summary(lm(barpm ~ barpred, data=spatialall))$r.squared)
res[res$type=="tempmean", 'm3.rmspe.space'] <- print(rmse(residuals(m1.fit.all.spat)))

#temporal
tempoall<-left_join(mod1,spatialall)
tempoall$delpm <-tempoall$temp.mean-tempoall$barpm
tempoall$delpred <-tempoall$pred.m3-tempoall$barpred
mod_temporal <- lm(delpm ~ delpred, data=tempoall)
res[res$type=="tempmean", 'm3.R2.time'] <-  print(summary(lm(delpm ~ delpred, data=tempoall))$r.squared)
saveRDS(res, "./data/outputs/AQUA/2002/c03/MEXICO.results.2002.AQ.tempmean.rds")



#create final prediction data set for use in health outcome studies

#import mod2.n
mod2.n<- readRDS( "./data/outputs/AQUA/2002/c03/MEXICO.mod2.2002.AQ.night.mean.predm2.rds")
mod2.n<-mod2.n[,c("lstid","day","pred.m2"),with=FALSE]

#----------------> store the best available
mod3best <- mod3[, list(lstid, long_lst, lat_lst, day, pred.m3)]
setkey(mod3best, day, lstid)
setkey(mod2.n, day, lstid)
mod3best <- merge(mod3best, mod2.n[,list(lstid, day, pred.m2)], all.x = T)
setkey(mod1,day,lstid)
mod3best <- merge(mod3best, mod1[,list(lstid,day,pred.m1,temp.mean)], all.x = T,allow.cartesian = T)
mod3best[,bestpred := pred.m3]
mod3best[!is.na(pred.m2),bestpred := pred.m2]
mod3best[!is.na(pred.m1),bestpred := pred.m1]
summary(mod3best$bestpred)
mod3best<-select(mod3best,day,lstid,long_lst,lat_lst,bestpred)
#save
saveRDS(mod3best,"./data/outputs/AQUA/2002/c03/MEXICO.2002.AQ.mean.bestpred.rds")
mod3best<-filter(mod3best,!is.na(bestpred))

#save for plotting in QGIS
out <- mod3best %>% group_by(lstid) %>%
  summarise(x=mean(long_lst, na.rm=TRUE), y =mean(lat_lst, na.rm=TRUE), bestpred=mean(bestpred, na.rm=TRUE))
out<-na.omit(out)
write.csv(out,"./data/outputs/AQUA/2002/c03/MEXICO.2002.AQ.mean.bestpredmap.csv")
#save res
saveRDS(res,"./data/outputs/AQUA/2002/c03/MEXICO.results.2002.AQ.tempmean.rds")

keep(rmse, sure=TRUE) 
gc()
