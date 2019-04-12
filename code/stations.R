suppressPackageStartupMessages(
   {library(car)
    library(data.table)
    library(dplyr)
    library(lubridate)
    library(plyr)
    library(readr)
    library(readxl)
    library(reshape2)
    library(stringi)
    library(purrr)})

data.root = "/data-belle/Mexico_temperature/stations/raw"

# Data from SIMAT
## REDMET

get.redmet = function()
   {l = sapply(c("TMP", "WSP", "RH", "pressure"), simplify = F, function(variable)
       {paths = list.files(
            paste0(data.root, (if (variable == "pressure")
                "/ATM_PRESSURE" else
                "/REDMET")),
            full.names = T,
            pattern = if (variable == "pressure")
               "*" else
               paste0("*", variable, ".xls"))
        d = rbindlist(lapply(paths, function(fname)
           {d = as.data.table(read_excel(fname))
            d[d == -99] <- NA
            setnames(d, c("FECHA", "HORA"), c("date", "time"))
            d = melt(d, id.vars = c("date","time"),
                variable.name = "stationid", value.name="value")
            d[,
               if (sum(!is.na(value)) >= 18)
                  {if (variable == "TMP")
                       .(
                            temp.mean = mean(value, na.rm = TRUE),
                            temp.max = max(value, na.rm = TRUE),
                            temp.min = min(value, na.rm = TRUE))
                   else
                        .(value.mean = mean(value, na.rm = TRUE))},
               by = .(date, stationid)]}))
        if (variable != "TMP" )
            setnames(d, "value.mean", (
                if (variable == "WSP") "wind.speed.mean"
                else if (variable == "RH") "relative.humidity.mean"
                else if (variable == "pressure") "pressure.mean"
                else stop()))
        d})
    merge(
        merge(
           merge(l$TMP, l$WSP, by = c("date", "stationid"), all = T),
           l$RH, by = c("date", "stationid"), all = T),
        l$pressure, by = c("date", "stationid"), all = T)}


## REDDA
get.redda.precipitation = function()
   {precipitation = as.data.table(rbindlist(lapply(
        list.files(paste0(data.root, "/REDDA"), full.names = T, pattern = "PPH.xls"),
        read_excel)))
    precipitation[precipitation == -99] <- NA
    setnames(precipitation, "FECHA", "date")
    melt(precipitation,
       id.vars = c("date"),
       variable.name = "stationid",
       value.name = "total.precipitation")}

# SMN

## OBSERVATORIES
get.smn = function()
   {smn_obs = fread(paste0(data.root, "/OBS_2018/OBS_HLY.CSV"))[
       `YEAR-MONTH-DAY` >= "2003-01-01"]
    colnames(smn_obs) = gsub("-", "_", colnames(smn_obs), fixed = T)
   
    setnames(smn_obs,
        c("Station_ID","ELEMENT_CODE","YEAR_MONTH_DAY"),
        c("stationid", "variable", "date"))
    smn_obs$variable <- car::recode(smn_obs$variable,
        "101 = 'temp'; 104 = 'pptn';105 = 'relhum'; 106 = 'atmpres'; else = 'other'")
    value_cols = grep("VALUE", colnames(smn_obs), val = T)
    smn_obs <- smn_obs[
        smn_obs$variable != "other",
        mget(c("stationid", "variable", "date", value_cols))]

    smn_obs[,
        (value_cols) := lapply(.SD, function(v)
           {v = as.numeric(v)
            v[v == -99999] = NA
            v}),
        .SDcols = value_cols]

    l = lapply(c("temp", "relhum", "atmpres", "pptn"), function(vname)
       {d <- melt(
            smn_obs[variable == vname,
                mget(c("stationid", "date", value_cols))],
            id.vars = c("date", "stationid"),
            variable.name = "time",
            value.name = "value")
        d = d[,
           if (sum(!is.na(value)) >= 18)
              {if (vname == "temp")
                   .(
                        temp.mean = mean(value, na.rm = TRUE),
                        temp.max = max(value, na.rm = TRUE),
                        temp.min = min(value, na.rm = TRUE))
               else if (vname == "pptn")
                    .(precipitation.total = sum(value, na.rm = TRUE))
               else
                    .(value.mean = mean(value, na.rm = TRUE))},
           by = .(date, stationid)]
        if (!vname %in% c("temp", "pptn"))
            setnames(d, "value.mean", (
                if (vname == "relhum") "relative.humidity.mean"
                else if (vname == "atmpres") "pressure.mean"
                else stop(vname)))
        d})

     Reduce(
       function(x, y)
           merge(x, y, all = T, by = c("stationid", "date")),
       l)}


#SMN
## ESIMES xls

smnesimesxls = list.files(path = paste0(data.root, "/ESIMEs_2018"),
   pattern="*.xlsx$", full.names = T)

smn_esimexls = lapply(smnesimesxls, function(fname){
   print(fname)
   d = read_excel(fname)
   if(any(colnames(d)=="EstaciÃ³n")){
       j = which(colnames(d)=="EstaciÃ³n")
       colnames(d)[j] = "Estacion"
   }
   d = read_excel(fname)
   if(any(colnames(d)=="Station")){
     j = which(colnames(d)=="Station")
     colnames(d)[j] = "Estacion"
   }
   d = read_excel(fname)
   if(any(colnames(d)=="DateTime")){
     j = which(colnames(d)=="DateTime")
     colnames(d)[j] = "Fecha-Tiempo"
   }
      if (all(colnames(d) %in% c("Estacion", "Fecha-Tiempo", "SR(W/m^2)", "Rain(mm) 10Min", "DP", "RH %", "GTempC(C) 12Hrs", "TS(C)", "ATCmax(C)", "ATCmin(C)", "ATC(C)", "WSmax(m/s)", "WSmin(m/s)", "WS(m/s)", "WDmax(m/s)", "WDmin(m/s)", "WD", "QFF", "BP(mbar)","VPS","VCC","Rain(mm) 24hrs","Rain(mm) 6Hrs","Rain(mm) 3Hrs","Rain(mm) 1Hr"))){
      data.table(d)
   } else {
          return(NULL)
  }
}) 

esimexls = rbindlist(smn_esimexls, use.names = TRUE, fill = TRUE)
setnames(esimexls, c("Estacion", "Fecha-Tiempo", "Rain(mm) 10Min", "RH %", "ATC(C)", "WS(m/s)", "BP(mbar)" ),
                   c("stationid", "date", "pptn", "relhum", "temp", "winds", "atmpres" ))
esimexls      <- esimexls %>% select("stationid", "date", "pptn", "relhum", "temp", "winds", "atmpres")
esimexls$date <- stri_replace(esimexls$date, fixed = "p.m.", replacement = "PM")
esimexls$date <- stri_replace(esimexls$date, fixed = "a.m.", replacement = "AM")
esimexls$date <- dmy_hms(esimexls$date)
esimexls$date <- with_tz(esimexls$date, tzone = "America/Mexico_City")
esimexls      <- melt(esimexls,
                      id.vars = c("date", "stationid"),
                      variable.name = "variable",
                      value.name="value")

esimexls[ , date := as.Date(date)]
esimexls[ , n := sum(!is.na(value)), by = .(date, stationid, variable)]
esimexls = na.omit(esimexls)
esimexls = esimexls[n>=108]

esimexls = list(
   esimexls[variable=="relhum"  ,.(relative.humidity   = mean(value)), by=.(date,stationid)]
  ,esimexls[variable=="winds"   ,.(wind.speed.mean     = mean(value)), by=.(date,stationid)]
  ,esimexls[variable=="atmpres" ,.(pressure.mean       = mean(value)), by=.(date,stationid)]
  ,esimexls[variable=="pptn"    ,.(precipitation.total = sum(value)),  by=.(date,stationid)]
  ,esimexls[variable=="temp"    ,.(temp.mean           = mean(value)                                                                                                            ,temp.min            = min(value)
                                  ,temp.max            = max(value)),  by=.(date,stationid)]
)
esimexls = reduce(esimexls,merge,all=TRUE)


## ESIMES csv

smnesimescsv = list.files(path = paste0(data.root, "/ESIMEs_2018"),
                          pattern="*.csv$", full.names = T)

smn_esimecsv = lapply(smnesimescsv, function(fname){
       names = read_csv(fname, n_max = 0) %>% names()
           d = read_csv(fname, col_names = names, skip = 2)
               setnames(d, c("Station", "DateTime", "Rh", "ATC", "Rain", "BP", "AvgWSV"),
                           c("stationid", "date", "relhum", "temp", "pptn", "atmpres", "winds"))
           d = melt(d, id.vars = c("date","stationid"),
           variable.name = "variable", value.name="value")
           d = transform(d, value = as.numeric(value),
                         variable = as.character(variable))
                variables <- c("relhum", "temp", "pptn", "atmpres", "winds")
                d[d$variable %in% variables,]
           d = as.data.table(d)
})

esimecsv = rbindlist(smn_esimecsv, use.names = TRUE, fill = TRUE)
esimecsv$date <- stri_replace(esimecsv$date, fixed = "p.m.", replacement = "PM")
esimecsv$date <- stri_replace(esimecsv$date, fixed = "a.m.", replacement = "AM")
esimecsv$date <- ymd_hms(esimecsv$date)
esimecsv$date <- with_tz(esimecsv$date, tzone = "America/Mexico_City")

esimecsv[ , date := as.Date(date)]
esimecsv[ , n := sum(!is.na(value)), by = .(date,stationid,variable)]
esimecsv = na.omit(esimecsv)
esimecsv = esimecsv[n>=108]

esimecsv = list(
                esimecsv[variable=="relhum"  ,.(relative.humidity   = mean(value)), by=.(date,stationid)]
               ,esimecsv[variable=="winds"   ,.(wind.speed.mean     = mean(value)), by=.(date,stationid)]
               ,esimecsv[variable=="atmpres" ,.(pressure.mean       = mean(value)), by=.(date,stationid)]
               ,esimecsv[variable=="pptn"    ,.(precipitation.total = sum(value)),  by=.(date,stationid)]
               ,esimecsv[variable=="temp"    ,.(temp.mean           = mean(value)                                                                                                        ,temp.min            = min(value)
                                               ,temp.max            = max(value)),  by=.(date,stationid)]
)

esimecsv = reduce(esimecsv,merge,all=TRUE)

smn_esimes <- as.data.table(rbind(esimexls, esimecsv))
rm(smn_esimexls,smn_esimecsv, esimexls, esimecsv)


##EMAS
smnemas = list.files(path = paste0(data.root, "/EMAs_2018"),
                     pattern="*.xls$", full.names = T)

smn_smnemas = lapply(smnemas, function(fname){
                      names = read_xls(fname, n_max = 0) %>% names()
                          d = read_xls(fname, col_names = names, skip = 2)
                       })
