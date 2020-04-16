####################
##
## Read data into dmort, dmorti, dhosp, dhospi:
##  dmort, dmorti: cumulative and incident deaths per day
##  dhosp, dhospi: cumulative and incident hospitalisations per day
##  dstartdate: first entry
##
####################

source("data.R")

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

##
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

death_underreporting_factor <- 1

dmort <- y.dmort + o.dmort
dmorti <- y.dmorti + o.dmorti

########################
## Age group demography 
########################

y.N <- 8.55E6
o.N <- N - y.N

## Group-specific factors for dead_rate with respect to estimated died_rate
y.died_rate_factor = 0.09412495/0.66
o.died_rate_factor = 2.34371716/0.66

## Group-specific factors for hosp_rate with respect to estimated hosp_rate
y.hosp_rate_factor = 0.09412495/0.66
o.hosp_rate_factor = 2.34371716/0.66
