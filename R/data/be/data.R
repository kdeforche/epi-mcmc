####################
##
## Read data into dmort, dmorti, dhosp, dhospi:
##  dmort, dmorti: cumulative and incident deaths per day
##  dhosp, dhospi: cumulative and incident hospitalisations per day
##  dstartdate: first entry
##
####################

dstartdate <- as.Date("2020/3/10")

be.hosp <- read.csv("COVID19BE_HOSP.csv")
print(aggregate(be.hosp$NEW_IN, by=list(date=be.hosp$DATE), FUN=sum))
dhospi <- aggregate(be.hosp$NEW_IN, by=list(date=be.hosp$DATE), FUN=sum)$x

## older data estimated from charts in older reports, missing in file
dhospi <- c(3, 12, 47, 50, 50, dhospi)

## remove latest data point, seems to be an incomplete day
## dhospi <- dhospi[1:(length(dhospi) - 1)]

dhosp <- cumsum(dhospi)

be.mort <- read.csv("COVID19BE_MORT.csv")
print(aggregate(be.mort$DEATHS, by=list(date=be.mort$DATE), FUN=sum))
dmorti <- aggregate(be.mort$DEATHS, by=list(date=be.mort$DATE), FUN=sum)$x

## remove two latest data point, may be incomplete because of WZC
dmorti <- dmorti[1:(length(dmorti) - 2)]
dmort <- cumsum(dmorti)

## all data series are now from 10/3

death_underreporting_factor <- 1

######
##
## Separate young and old deaths
##
######

be.mort$y = (be.mort$AGEGROUP == "0-24" | be.mort$AGEGROUP == "25-44" | be.mort$AGEGROUP == "45-64")
be.mort.na = subset(be.mort, is.na(be.mort$y))
na.dmorti = aggregate(be.mort.na$DEATHS, by=list(date=be.mort.na$DATE), FUN=sum, drop=F)
na.dmorti$x[is.na(na.dmorti$x)] = 0

ycount = sum(be.mort$y == T & as.character(be.mort$DATE) > "2020-03-20", na.rm=TRUE)
ocount = sum(be.mort$y == F & as.character(be.mort$DATE) > "2020-03-20", na.rm=TRUE)
yfract = ycount / (ycount + ocount)
yfract = yfract * 0.2 # assume that maybe in particular older people do not have their age properly

mortaggr = aggregate(be.mort$DEATHS, by=list(y=be.mort$y, date=be.mort$DATE), FUN=sum, drop=F)
mortaggr$x[is.na(mortaggr$x)] = 0

y.dmorti = mortaggr$x[mortaggr$y == T]
o.dmorti = mortaggr$x[mortaggr$y == F]

y.dmorti = y.dmorti + floor(yfract * na.dmorti$x)
o.dmorti = o.dmorti + floor((1 - yfract) * na.dmorti$x)

## estimated from the epidemiological report
## 

o.wzcmorti = c(0, 0, 0, 0, 2, 1, 1, 2, 4, 10, 9, 12, 21, 19, 12, 24, 25,
               40, 44, 45, 48, 53, 100, 70, 103, 127, 102, 146, 179, 205,
               160, 192, 150, 180, 10, 2)

print(paste("last day o.wzcmorti: ", dstartdate + length(o.wzcmorti) - 1))

o.dmorti = o.dmorti - o.wzcmorti

## remove two latest data point, may be incomplete because of WZC
y.dmorti <- y.dmorti[1:(length(y.dmorti) - 2)]
o.dmorti <- o.dmorti[1:(length(o.dmorti) - 2)]

y.dmort <- cumsum(y.dmorti)
o.dmort <- cumsum(o.dmorti)

par(nrow=3)
barplot(y.dmorti + o.dmorti)
barplot(y.dmorti)
barplot(o.dmorti)

print(paste("last day morti: ", dstartdate + length(o.dmorti) - 1))

#####################
## Population size
#####################

N <- 11.5E6

#####################
## Lockdown measures
#####################

## Date of lockdown phase (1)
lockdown_offset <- as.numeric(as.Date("2020/3/13") - dstartdate)

## over how many days the lockdown is estimated to have occurred
lockdown_transition_period <- 7

## how many deaths at date of lockdown
total_deaths_at_lockdown <- 5
