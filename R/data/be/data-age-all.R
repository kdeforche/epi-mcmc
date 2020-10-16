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

y.dmorti <- y.dmorti[1:(length(y.dmorti) - 2)]
o.dmorti <- o.dmorti[1:(length(o.dmorti) - 2)]

y.dmort <- cumsum(y.dmorti)
o.dmort <- cumsum(o.dmorti)

par(mfrow=c(1,3))
barplot(y.dmorti)
barplot(o.dmorti)
barplot(y.dmorti + o.dmorti)

o.wdmorti <- rep(1, length(o.dmorti))
o.wdmorti[154:163] = 0.1

print(paste("last day morti: ", dstartdate + length(o.dmorti) - 1))

dmort <- y.dmort + o.dmort
dmorti <- y.dmorti + o.dmorti



##
## time series of IFR estimates per group

be.mort$week = format(as.Date(be.mort$DATE), "%Y-%V")

lastweek <- be.mort$week[length(be.mort$week)]
if (length(be.mort$week[be.mort$week == lastweek]) < 10) {
    be.mort = subset(be.mort, week != lastweek)
}

weekgroup = aggregate(be.mort$DEATHS, by=list(week=be.mort$week,group=be.mort$AGEGROUP), FUN=sum, drop=F)
weekgroup$x[is.na(weekgroup$x)] = 0


##g.ifr <- c(0.009604, 0.090175, 0.82025, 3.105, 6.04, 11.7) / 100
## Based on "Belgian Covid-19 Mortaility ...", Geert Molenbergs et. al; Table 6
g.ifr <- c(0.0005, 0.017, 0.21, 2.2, 4.29, 11.8) / 100

ifrf = 0
ifrd = as.Date("2020-09-15")

calcifr.y <- function(x) {
    x <- x + ifr_epsilon ## As if a multinomial prior, would be better if we do some windowing

    if (i > (ifrd - dstartdate)) {
        x <- x + ifrf * ifr_epsilon
    }
    i <<- i + 7
    
    y.ifr <- sum(x[1:3]) / sum(x[1:3] / g.ifr[1:3])
    y.ifr
}   

calcifr.o <- function(x) {
    x <- x + ifr_epsilon ## As if a multinomial prior

    o.ifr <- sum(x[4:6]) / sum(x[4:6] / g.ifr[4:6])
    o.ifr
}   

calcifr <- function(x) {
    x <- x + ifr_epsilon  ## As if a multinomial prior

    if (i > (ifrd - dstartdate)) {
        x[1:3] <- x[1:3] + ifrf * ifr_epsilon
    }
    i <<- i + 7

    ifr <- sum(x) / sum(x / g.ifr)
    ifr
}   

i <<- 1
y.weekifr = aggregate(weekgroup$x, by=list(week=weekgroup$week), FUN=calcifr.y, drop=F)
o.weekifr = aggregate(weekgroup$x, by=list(week=weekgroup$week), FUN=calcifr.o, drop=F)
i <<- 1
all.weekifr = aggregate(weekgroup$x, by=list(week=weekgroup$week), FUN=calcifr, drop=F)

weekifr <- data.frame(y.weekifr$week, y.weekifr$x, o.weekifr$x, all.weekifr$x)
weekifr$index = (seq(3,length(o.dmorti) + 7,7))[1:(length(y.weekifr$week))]

m.yifr <- smooth.spline(weekifr$index, weekifr$y.weekifr.x, df=5)
yf <- data.frame(index=seq(1:length(o.dmorti)))
y.ifr <- as.numeric(unlist(predict(m.yifr, yf)$y))

m.oifr <- smooth.spline(weekifr$index, weekifr$o.weekifr.x, df=3)
o.ifr <- as.numeric(unlist(predict(m.oifr, yf)$y))

m.ifr <- smooth.spline(weekifr$index, weekifr$all.weekifr.x, df=6)
all.ifr <- as.numeric(unlist(predict(m.ifr, yf)$y))

pdf("ifr.pdf", width=15, height=5)
par(mfrow=c(1,3))
x <-seq(dstartdate, dstartdate+length(y.ifr)-1, by=1)
plot(x, y.ifr * 100, type='l', main="Time profile of COVID-19 IFR Belgium (<65y)",
        xlab="Date", ylab="Infection Fatality Rate (%)", ylim=c(0, 0.1))
points(dstartdate + weekifr$index, weekifr$y.weekifr.x * 100)
lines(x, y.ifr * 100, type='l', col=2)

plot(x, o.ifr * 100, type='l', main="Time profile of COVID-19 IFR Belgium (>65y)",
     xlab="Date", ylab="Infection Fatality Rate (%)", ylim=c(0, 7))
points(dstartdate + weekifr$index, weekifr$o.weekifr.x * 100)

plot(x, all.ifr * 100, type='l', main="Time profile of COVID-19 IFR Belgium (log scale)",
     xlab="Date", ylab="Infection Fatality Rate (%)", ylim=c(0.02, 2), log="y")
points(dstartdate + weekifr$index, weekifr$all.weekifr.x * 100)

print(paste(c("Current IFR (%): ", all.ifr[length(all.ifr)-1] * 100)))

dev.off()

### HR
### Statistics for Belgium  up to June 14
### https://covid-19.sciensano.be/sites/default/files/Covid19/COVID-19_THEMATIC%20REPORT_COVID-19%20HOSPITALISED%20PATIENTS_NL.pdf
###
hosppct.10 <- c(1.2, 0.5, 2.1, 4.1, 8.0, 15.0, 17.0, 20.7, 24.1, 7.3)
hosppct.6 <- c(1.2 + 0.5 + 2.1/2, 2.1/2 + 4.1 + 8/2, 8/2 + 15 + 17/2, 17/2 + 20.7/2, 20.7/2 + 24.1/2, 24.1/2 + 7.3) * 15139 / 100

mortaggr.2 = aggregate(be.mort$DEATHS, by=list(group=be.mort$AGEGROUP, date=be.mort$DATE), FUN=sum, drop=FALSE)
mortaggr.2$x[is.na(mortaggr.2$x)] = 0
mortaggr.3 <- subset(mortaggr.2, as.Date(date) < as.Date("2020-06-14"))

mort.june14 <- aggregate(mortaggr.3$x, by=list(group=mortaggr.3$group), FUN=sum, drop=FALSE)

g.hr <- hosppct.6 / (mort.june14$x / g.ifr)
g.hr / g.ifr

require(ggplot2)
require(reshape2)

rates <- data.frame(group = levels(be.mort$AGEGROUP), hr = g.hr * 100, ifr = g.ifr * 100)
rates.1 <- melt(rates, id.vars = c("group"), value.name="value")
rates.1$variable = as.character(rates.1$variable)
rates.1$variable[rates.1$variable=="hr"] = "Hospitalization rate (HR)"
rates.1$variable[rates.1$variable=="ifr"] = "Infection fatality rate (IFR)"
pdf("hrifr.pdf", width=6, height=3)
ggplot(data=rates.1, aes(x=group, y=value, fill=variable)) + geom_bar(stat="identity", position=position_dodge()) + scale_fill_discrete(name = "") + xlab("Age group") + ylab("Rate (%)")
dev.off()

calchr.y <- function(x) {
    x <- x + ifr_epsilon ## As if a multinomial prior, would be better if we do some windowing

    if (i > (ifrd - dstartdate)) {
         x <- x + ifrf * ifr_epsilon
    }
    i <<- i + 7

    inf <- x[1:3] / g.ifr[1:3]
    hosp <- g.hr[1:3] * inf
    
    y.hr <- sum(hosp) / sum(inf)

    y.hr
}   

calchr.o <- function(x) {
    x <- x + ifr_epsilon ## As if a multinomial prior

    inf <- x[4:6] / g.ifr[4:6]
    hosp <- g.hr[4:6] * inf

    o.hr <- sum(hosp) / sum(inf)
    o.hr
}   

calchr <- function(x) {
    x <- x + ifr_epsilon  ## As if a multinomial prior

    if (i > (ifrd - dstartdate)) {
         x[1:3] <- x[1:3] + ifrf * ifr_epsilon
    }
    i <<- i + 7
    
    inf <- x / g.ifr
    hosp <- g.hr * inf

    hr <- sum(hosp) / sum(inf)
    hr
}   

i <<- 1
y.weekhr = aggregate(weekgroup$x, by=list(week=weekgroup$week), FUN=calchr.y, drop=F)
o.weekhr = aggregate(weekgroup$x, by=list(week=weekgroup$week), FUN=calchr.o, drop=F)
i <<- 1
all.weekhr = aggregate(weekgroup$x, by=list(week=weekgroup$week), FUN=calchr, drop=F)

weekhr <- data.frame(y.weekhr$week, y.weekhr$x, o.weekhr$x, all.weekhr$x)
weekhr$index = (seq(3,length(o.dmorti) + 7,7))[1:(length(y.weekhr$week))]

m.yhr <- smooth.spline(weekhr$index, weekhr$y.weekhr.x, df=5)
yf <- data.frame(index=seq(1:length(o.dmorti)))
y.hr <- as.numeric(unlist(predict(m.yhr, yf)$y))

m.ohr <- smooth.spline(weekhr$index, weekhr$o.weekhr.x, df=3)
o.hr <- as.numeric(unlist(predict(m.ohr, yf)$y))

m.hr <- smooth.spline(weekhr$index, weekhr$all.weekhr.x, df=4)
all.hr <- as.numeric(unlist(predict(m.hr, yf)$y))

pdf("hr.pdf", width=15, height=5)
par(mfrow=c(1,3))
x <-seq(dstartdate, dstartdate+length(y.hr)-1, by=1)
plot(x, y.hr * 100, type='l', main="Time profile of COVID-19 HR Belgium (<65y)",
        xlab="Date", ylab="Hospitalization Rate (%)", ylim=c(0, 1))
points(dstartdate + weekhr$index, weekhr$y.weekhr.x * 100)

plot(x, o.hr * 100, type='l', main="Time profile of COVID-19 HR Belgium (>65y)",
     xlab="Date", ylab="Hospitalization Rate (%)", ylim=c(0, 10))
points(dstartdate + weekhr$index, weekhr$o.weekhr.x * 100)

plot(x, all.hr * 100, type='l', main="Time profile of COVID-19 HR Belgium (log scale)",
     xlab="Date", ylab="Hospitalization Rate (%)", log="y", ylim=c(0.1, 10))
points(dstartdate + weekhr$index, weekhr$all.weekhr.x * 100)

print(paste(c("Current HR (%): ", all.hr[length(all.hr)-1] * 100)))

dev.off()


#####
## Estimated number of total infected per age group from IFR:
#####

g.N <- c(3228894, 2956684, 3080528, 1147009, 690685, 327606)
g.aggr <- aggregate(be.mort$DEATHS, by=list(group=be.mort$AGEGROUP), FUN=sum, drop=F)
g.deaths <- g.aggr$x
print(g.deaths)
g.infected <- g.deaths / g.ifr
g.pct.infected <- g.infected / g.N * 100

barplot(g.pct.infected ~ g.aggr$group)
print(sum(g.infected[1:3]) / sum(g.N[1:3]))
print(sum(g.infected[4:6]) / sum(g.N[4:6]))

########################
## Age group demography 
########################

y.N <- 8.55E6
o.N <- N - y.N

###################
##
###################

glogis <- function(t, A, K, C, Q, B, M, v) {
    A + (K - A)/((C + Q * exp(-B * (t - M)))^(1/v))
}

vs <- 1:length(y.ifr)

g <- glogis(vs, 0, 1, 1, 0.5, 0.1, 200, 0.5)
##g <- glogis(vs, 0, 1, 1, 0.5, 0.05, 140, 0.5)

plot(x, y.ifr, type='l', ylim=c(0, 1E-3))
points(x, y.ifr * (1 - 0.5 * g))

##g1 <- glogis(vs, 0, 1, 1, 0.5, 0.05, 70, 0.5)
##plot(x, o.hr, type='l', ylim=c(0, 0.3))
## o.hr <- o.hr * (1 + 3 * g1)
