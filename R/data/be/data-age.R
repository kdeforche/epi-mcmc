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
## hospi : can be approximately done using statistics from the epi report
##

## date 0-49 50-79 80+
## 08-14/3 : 20 55 25
## 15-21/3 : 18 58 25
## 22-28/3 : 17 58 26
## 29-04/4 : 15 56 29
## 05-11/4 : 13 51 36
## 12-18/4 : 13 46 40
##
## 50-59: 15.4
## 60-69: 18.0
## 70-79: 21.0
##  < 65 : 24.4; > 65 : 30
##   compare y.dhospi en o.dhospi apart
## -> should increase hospitalisations for < 65 vs for > 65y and thus
##    predict higher number of 

hosp_week_age_statistics <- data.frame(r0.49 = c(20, 18, 17, 15, 14, 13, 18, 18, 14))
hosp_week_age_statistics$r50.79 = c(55, 59, 59, 56, 51, 47, 46, 47, 43)
hosp_week_age_statistics$r80 = c(25, 25, 26, 29, 36, 41, 38, 37, 45)

for (i in 1:dim(hosp_week_age_statistics)[1]) {
    hosp_week_age_statistics[i,] <- hosp_week_age_statistics[i,] / sum(hosp_week_age_statistics[i,])
}

bucket_50.79.y = 24.4/(24.4 + 30)
hosp_week_age_statistics$y = hosp_week_age_statistics[,1] + hosp_week_age_statistics[,2] * bucket_50.79.y
hosp_week_age_statistics$index = 1:dim(hosp_week_age_statistics)[1] * 7 - 5

model <- lm(y ~ poly(index, 3), data=hosp_week_age_statistics)
yf <- data.frame(index=seq(1:length(dhospi)))

pdf("hosp-age.pdf", width=8, height=6)

plot(dstartdate + (yf$index - 1), predict(model, yf) * 100, type='l', ylim=c(0, 100), xlab="Date", ylab="Fraction of hospitalisations", main="Distribution of hospitalisations in younger and older groups", col='blue')
points(dstartdate + (hosp_week_age_statistics$index - 1), hosp_week_age_statistics$y * 100, col='blue')
lines(dstartdate + (yf$index - 1), 100 - predict(model, yf) * 100, type='l', col='darkgreen')
points(dstartdate + (hosp_week_age_statistics$index - 1), 100 - hosp_week_age_statistics$y * 100, col='darkgreen')
legend("topleft", inset=0.02, legend=c("< 65", ">= 65"),
       col=c("blue", "darkgreen"),lty=1)

dev.off()

y.dhospi = round(dhospi * predict(model, yf))
o.dhospi = dhospi - y.dhospi

##
## dmorti : ages are recorded, substracting wzc needs estimates from the epi report

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

o.wzcmorti = c(0, 0, 0, 0, 2, 1, 1, 2, 3, 10, 9, 12, 21, 19, 12, 24, 25,
               40, 49, 52, 48, 60, 122, 76, 112, 137, 110, 156, 182, 185,
               170, 199, 176, 240, 155, 165, 149, 153, 117, 130, 126, 120, # 20/04/20
               116, 112, 95, 74, 76, 95, 70, 55, 48, 46, 47, 33, 45, 40, 38, # 05/05/20
               37, 40, 28, 35, 38, 25, 24, 20, 23, 11, 8, 21, 12, 0, 0) # 20/05/20

print(paste("last day o.wzcmorti: ", dstartdate + length(o.wzcmorti) - 1))
print(paste("last day o.dmorti: ", dstartdate + length(o.dmorti) - 1))

o.dmorti = o.dmorti - o.wzcmorti

## remove two latest data point, may be incomplete because of WZC
y.dmorti <- y.dmorti[1:(length(y.dmorti) - 2)]
o.dmorti <- o.dmorti[1:(length(o.dmorti) - 2)]

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
