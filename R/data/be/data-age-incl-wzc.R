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

##
## dcasei : ages are recorded

levs <- as.character(seq(min(as.Date(be.case$DATE[!is.na(be.case$DATE)])), max(as.Date(be.case$DATE[!is.na(be.case$DATE)])), by=1))
dates <- data.frame(date=levs)

be.case$y = (be.case$AGEGROUP == "0-9" | be.case$AGEGROUP == "10-19" | be.case$AGEGROUP == "20-29" | be.case$AGEGROUP == "30-39" | be.case$AGEGROUP == "40-49" | be.case$AGEGROUP == "50-59")
be.case.na = subset(be.case, is.na(be.case$y))
na.dcasei.1 = aggregate(be.case.na$CASES, by=list(date=be.case.na$DATE), FUN=sum, drop=F)
na.dcasei = merge(dates, na.dcasei.1, by.x=c("date"), by.y=c("date"), all=TRUE)
na.dcasei$x[is.na(na.dcasei$x)] = 0

ycount = sum(be.case$y == T, na.rm=TRUE)
ocount = sum(be.case$y == F, na.rm=TRUE)
yfract = ycount / (ycount + ocount)

caseaggr.1 = aggregate(be.case$CASES, by=list(y=be.case$y, date=be.case$DATE), FUN=sum, drop=F)
caseaggr = merge(dates, caseaggr.1, by.x=c("date"), by.y=c("date"), all=TRUE)
caseaggr$x[is.na(caseaggr$x)] = 0

y.dcasei = caseaggr$x[caseaggr$y == T]
o.dcasei = caseaggr$x[caseaggr$y == F]

y.dcasei = y.dcasei + floor(yfract * na.dcasei$x)
o.dcasei = o.dcasei + floor((1 - yfract) * na.dcasei$x)

y.dcasei <- y.dcasei[10:length(y.dcasei)]
o.dcasei <- o.dcasei[10:length(o.dcasei)]

y.dcasei <- y.dcasei[1:(length(y.dcasei) - 2)]
o.dcasei <- o.dcasei[1:(length(o.dcasei) - 2)]

y.dcase <- cumsum(y.dcasei)
o.dcase <- cumsum(o.dcasei)

dcase <- y.dcase + o.dcase
dcasei <- y.dcasei + o.dcasei

y.dhospi = y.dcasei
o.dhospi = o.dcasei
y.dhosp = y.dcase
o.dhosp = o.dcase
dhospi = dcasei
dhosp = dcase

##
## dmorti : ages are recorded

levs <- as.character(seq(min(as.Date(be.mort$DATE)), max(as.Date(be.mort$DATE)), by=1))
dates <- data.frame(date=levs)
df.2 <- transform(be.mort, DATE=factor(DATE, sort(unique(levs))))
be.mort <- df.2
be.mort$y = (be.mort$AGEGROUP == "0-24" | be.mort$AGEGROUP == "25-44" | be.mort$AGEGROUP == "45-64")
be.mort.na = subset(be.mort, is.na(be.mort$y))
na.dmorti.1 = aggregate(be.mort.na$DEATHS, by=list(date=be.mort.na$DATE), FUN=sum, drop=FALSE)
na.dmorti = merge(dates, na.dmorti.1, by.x=c("date"), by.y=c("date"), all=TRUE)
na.dmorti$x[is.na(na.dmorti$x)] = 0

ycount = sum(be.mort$y == T & as.character(be.mort$DATE) > "2020-03-20", na.rm=TRUE)
ocount = sum(be.mort$y == F & as.character(be.mort$DATE) > "2020-03-20", na.rm=TRUE)
yfract = ycount / (ycount + ocount)
yfract = yfract * 0.2 # assume that maybe in particular older people do not have their age properly

mortaggr.1 = aggregate(be.mort$DEATHS, by=list(y=be.mort$y, date=be.mort$DATE), FUN=sum, drop=FALSE)
mortaggr = merge(dates, mortaggr.1, by.x=c("date"), by.y=c("date"), all=TRUE)
mortaggr$x[is.na(mortaggr$x)] = 0

y.dmorti = mortaggr$x[mortaggr$y == T]
o.dmorti = mortaggr$x[mortaggr$y == F]

y.dmorti[is.na(y.dmorti)] = 0
o.dmorti[is.na(o.dmorti)] = 0

y.dmorti = y.dmorti + floor(yfract * na.dmorti$x)
o.dmorti = o.dmorti + floor((1 - yfract) * na.dmorti$x)

##y.dmorti <- y.dmorti[1:(length(y.dmorti)]
##o.dmorti <- o.dmorti[1:(length(o.dmorti)]

y.dmort <- cumsum(y.dmorti)
o.dmort <- cumsum(o.dmorti)

par(nrow=3)
barplot(y.dmorti)
barplot(o.dmorti)
barplot(y.dmorti + o.dmorti)

print(paste("last day morti: ", dstartdate + length(o.dmorti) - 1))

dmort <- y.dmort + o.dmort
dmorti <- y.dmorti + o.dmorti

##
## time series of IFR estimates per group

be.mort$week = format(as.Date(be.mort$DATE), "%Y-%V")
weekgroup = aggregate(be.mort$DEATHS, by=list(week=be.mort$week,group=be.mort$AGEGROUP), FUN=sum, drop=F)
weekgroup$x[is.na(weekgroup$x)] = 0

g.ifr <- c(0.009604, 0.090175, 0.82025, 3.105, 6.04, 11.7) / 100

calcifr.y <- function(x) {
    x <- x + 1 ## As if a multinomial prior

    y.ifr <- sum(x[1:3]) / sum(x[1:3] / g.ifr[1:3])
    y.ifr
}   

calcifr.o <- function(x) {
    x <- x + 1 ## As if a multinomial prior

    o.ifr <- sum(x[4:6]) / sum(x[4:6] / g.ifr[4:6])
    o.ifr
}   

calcifr <- function(x) {
    x <- x + 1 ## As if a multinomial prior

    ifr <- sum(x) / sum(x / g.ifr)
    ifr
}   

y.weekifr = aggregate(weekgroup$x, by=list(week=weekgroup$week), FUN=calcifr.y, drop=F)
o.weekifr = aggregate(weekgroup$x, by=list(week=weekgroup$week), FUN=calcifr.o, drop=F)
all.weekifr = aggregate(weekgroup$x, by=list(week=weekgroup$week), FUN=calcifr, drop=F)

weekifr <- data.frame(y.weekifr$week, y.weekifr$x, o.weekifr$x, all.weekifr$x)
weekifr$index = seq(3,length(o.dmorti),7)

m.yifr <- smooth.spline(weekifr$index, weekifr$y.weekifr.x, df=6)
yf <- data.frame(index=seq(1:length(o.dmorti)))
y.ifr <- as.numeric(unlist(predict(m.yifr, yf)$y))

m.oifr <- smooth.spline(weekifr$index, weekifr$o.weekifr.x, df=6)
o.ifr <- as.numeric(unlist(predict(m.oifr, yf)$y))

m.ifr <- smooth.spline(weekifr$index, weekifr$all.weekifr.x, df=6)
all.ifr <- as.numeric(unlist(predict(m.ifr, yf)$y))

plot(y.ifr, type='l')
points(weekifr$index, weekifr$y.weekifr.x)

plot(o.ifr, type='l')
points(weekifr$index, weekifr$o.weekifr.x)

########################
## Age group demography 
########################

y.N <- 8.55E6
o.N <- N - y.N
