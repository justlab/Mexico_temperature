library(magrittr)
library(plyr)
library(reshape2)
library(ggplot2)

#################################################################################
#################################################################################
#################################################################################
# Prepare data

# setwd("/media/qnap_eea/Projects/P018.IL.Israel.LST.Ta/2.work")
setwd("R:\\work")

files = list.files(pattern = "\\.CV\\.rds")

result = NULL

for(i in files) {

  dat = readRDS(i)
  
  s = strsplit(i, "\\.")
  s = s[[1]]
  
  if(s[6] == "max") dat$obs = dat$predmax
  if(s[6] == "mean") dat$obs = dat$predmean
  if(s[6] == "min") dat$obs = dat$predmin

  dat$sat = s[4]
  dat$var = s[6]
  
  dat = select(dat,stn, lat_lst, long_lst, day, sat, var, obs, pred.m1.cv)
  #dat = dat[, c("stn", "lat_lst", "long_lst", "day", "sat", "var", "obs", "pred.m1.cv")]
  dat <-as.data.table(dat)
  result = rbind(result, dat , fill=TRUE)
  
}

result$resid = result$obs - result$pred.m1.cv

#################################################################################
#################################################################################
#################################################################################
# Plots

# setwd("/home/michael/Dropbox/BGU/Adar/paper_graphs")

#################################################################################
# Distribution of residuals in sat / var

dat = result

ggplot(dat, aes(x = var, y = resid, fill = sat)) + 
  geom_boxplot() 

# Map

dat = result
dat$resid = abs(dat$resid)

dat = 
  ddply(
    dat,
    c("stn", "sat", "var"),
    summarize,
    lon = mean(long_lst),
    lat = mean(lat_lst),
    resid_min = min(resid),
    resid_max = max(resid),
    resid_50 = quantile(resid, 0.5)
  )
dat$type = paste0(dat$sat, "-", dat$var)
dat$type = factor(dat$type, levels = c("AQ-min", "AQ-mean", "AQ-max", "TR-min", "TR-mean", "TR-max"))
dat = melt(dat, measure.vars = c("resid_min", "resid_50", "resid_max"))

ggplot(dat, aes(x = lon, y = lat)) +
  geom_point(shape = 1, aes(size = value, colour = variable)) +
  scale_size_area(max_size = 20, breaks = c(0.1, 0.5, 1, 5, 10)) +
  scale_color_manual(values = c("blue", "black", "red")) +
  facet_wrap(~ type, nrow = 1) +
  coord_map()

# To add: separate maps for min/max/mean with color specifying sign of residual

#################################################################################
# Time-Series

dat = result
dat$day = paste0("2015", substr(dat$day, 5, 10))
dat$day = as.Date(dat$day)

dat = 
  ddply(
    dat,
    c("day", "sat", "var"),
    summarize,
    resid_min = min(resid),
    resid_max = max(resid),
    resid_05 = quantile(resid, 0.05),
    resid_95 = quantile(resid, 0.95)
  )
dat$type = paste0(dat$sat, "-", dat$var)
dat$type = factor(dat$type, levels = c("AQ-min", "AQ-mean", "AQ-max", "TR-min", "TR-mean", "TR-max"))
dat = melt(dat, measure.vars = c("resid_min", "resid_05", "resid_95", "resid_max"))

Sys.setlocale(locale = "C")

ggplot(dat, aes(x = day, y = value)) +
  geom_line(aes(colour = variable)) +
  scale_color_manual(values = c("blue", "black", "black", "red")) +
  scale_x_date(date_labels = "%b") +
  facet_wrap(~ type, ncol = 1)

dat$group = paste(dat$type, dat$variable)
ggplot(dat, aes(x = day, y = value)) +
  geom_line(aes(colour = var, group = group)) +
  # scale_color_manual(values = c("blue", "black", "black", "red")) +
  scale_x_date(date_labels = "%b") +
  facet_wrap(~ sat, ncol = 1)


#################################################################################
# Time-Series - Monthly

dat = result
dat$month = 
  dat$day %>% 
  as.Date %>% 
  format(format = "%b") %>% 
  factor(levels = month.abb)

ggplot(dat, aes(x = month, y = resid, colour = var)) +
  geom_boxplot() +
  # coord_cartesian(ylim = c(-1.5, 1.5)) +
  facet_grid(sat ~ .)

#################################################################################

View(result)
















