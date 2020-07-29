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

be.case$y = (be.case$AGEGROUP == "0-9" | be.case$AGEGROUP == "10-19" | be.case$AGEGROUP == "20-29" | be.case$AGEGROUP == "30-39" | be.case$AGEGROUP == "40-49" | be.case$AGEGROUP == "50-59")
be.case.na = subset(be.case, is.na(be.case$y))
na.dcasei = aggregate(be.case.na$CASES, by=list(date=be.case.na$DATE), FUN=sum, drop=F)
na.dcasei$x[is.na(na.dcasei$x)] = 0

ycount = sum(be.case$y == T, na.rm=TRUE)
ocount = sum(be.case$y == F, na.rm=TRUE)
yfract = ycount / (ycount + ocount)

caseaggr = aggregate(be.case$CASES, by=list(y=be.case$y, date=be.case$DATE), FUN=sum, drop=F)
caseaggr$x[is.na(caseaggr$x)] = 0

y.dcasei = caseaggr$x[caseaggr$y == T]
o.dcasei = caseaggr$x[caseaggr$y == F]

y.dcasei = y.dcasei + floor(yfract * na.dcasei$x)
o.dcasei = o.dcasei + floor((1 - yfract) * na.dcasei$x)

y.dcasei <- y.dcasei[10:length(y.dcasei)]
o.dcasei <- o.dcasei[10:length(o.dcasei)]

y.dcasei <- y.dcasei[1:(length(y.dcasei) - 1)]
o.dcasei <- o.dcasei[1:(length(o.dcasei) - 1)]

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
df.2 <- transform(be.mort, DATE=factor(DATE, sort(unique(levs))))
be.mort <- df.2
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

########################
## Age group demography 
########################

y.N <- 8.55E6
o.N <- N - y.N
