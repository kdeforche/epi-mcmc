print(c("Country:", country3))

ined <- read.csv("covid_pooled_01_09.csv")

country.data <- subset(ined, country_code==country3)

agegroups <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90+")
s <- subset(country.data, age_group %in% agegroups & excelsheet=="PHAS_Data")
s$date = as.Date(s$death_reference_date, format="%d.%m.%Y")

dstartdate <- min(s$date)

## We need to fill in the first data points
dmort <- c(1, 1, 1, 1, 3, 7, 8, 10, 10, 16, 20, 21, 25, 36, 42, 66, 92, 102, 110, 146, 180)
fract <- head(s, 10)$cum_death_both / sum(head(s, 10)$cum_death_both)

Round <- function(x, target) {
  r.x <- round(x)
  diff.x <- round(x) - x
  if ((s <- sum(r.x)) == target) {
    return(r.x)
  } else if (s > target) {
    select <- seq(along=x)[diff.x > 0]
    which <- which.max(diff.x[select])
    x[select[which]] <- r.x[select[which]] - 1
    Round(x, target)
  } else {
    select <- seq(along=x)[diff.x < 0]
    which <- which.min(diff.x[select])
    x[select[which]] <- r.x[select[which]] + 1
    Round(x, target)
  }
}

d <- c()
for (m in dmort) {
    d <- c(d, Round(m * fract, m))
}

s0dates <- seq(dstartdate - length(dmort), dstartdate - 1, 1)
s0 <- data.frame(pop_both = rep(1E20, length(d)), date = rep(s0dates, each=10), age_group = rep(agegroups, length(dmort)), cum_death_both = d)

s1 <- subset(s, select=c("pop_both", "date", "age_group", "cum_death_both"))

s = rbind(s0, s1)

levs <- seq(min(s$date), max(s$date), by=1)
dates <- data.frame(date=levs)
dstartdate <- min(s$date)

s$y = s$age_group %in% c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59")

mortaggr.1 = aggregate(s$cum_death_both, by=list(y=s$y, date=s$date), FUN=sum, drop=FALSE)
mortaggr <- merge(dates, mortaggr.1, by.x=c("date"), by.y=c("date"), all=TRUE)
mortaggr$x[is.na(mortaggr$x)] = 0

y.dmort = mortaggr$x[mortaggr$y == T]
o.dmort = mortaggr$x[mortaggr$y == F]

for (i in 1:length(y.dmort)) {
    if (is.na(y.dmort[i]))
        y.dmort[i] = y.dmort[i-1]
    if (is.na(o.dmort[i]))
        o.dmort[i] = o.dmort[i-1]
}

y.dmorti = pmax(diff(y.dmort), 0)
o.dmorti = pmax(diff(o.dmort), 0)

dmorti = y.dmorti + o.dmorti
dmort = y.dmort + o.dmort

s$death = s$cum_death_both

for (g in agegroups) {
    s$death[s$age_group == g] = c(s$cum_death_both[s$age_group == g][1],
                                  pmax(0, diff(s$cum_death_both[s$age_group == g])))
}

s$week = format(s$date, "%Y-%V")
s$group = as.factor(as.character(s$age_group))

weekgroup = aggregate(s$death, by=list(week=s$week,group=s$group), FUN=sum, drop=F)

## https://www.medrxiv.org/content/10.1101/2020.08.24.20180851v1 table S4
## https://www.folkhalsomyndigheten.se/contentassets/53c0dc391be54f5d959ead9131edb771/infection-fatality-rate-covid-19-stockholm-technical-report.pdf
g.ifr <- c(0.0005, 0.001, (0.004 + 0.009)/2, (0.017 + 0.029)/2, (0.053 + 0.086)/2,
           (0.154 + 0.241)/2, (0.359 + 0.642)/2, 1.92, 7.2, 16.2) / 100

calcifr.y <- function(x) {
    x <- x + ifr_epsilon ## As if a multinomial prior, would be better if we do some windowing

    y.ifr <- sum(x[1:6]) / sum(x[1:6] / g.ifr[1:6])
    y.ifr
}   

calcifr.o <- function(x) {
    x <- x + ifr_epsilon ## As if a multinomial prior, would be better if we do some windowing

    y.ifr <- sum(x[7:9]) / sum(x[7:9] / g.ifr[7:9])
    y.ifr
}   

calcifr <- function(x) {
    x <- x + ifr_epsilon ## As if a multinomial prior, would be better if we do some windowing

    ifr <- sum(x) / sum(x / g.ifr)
    ifr
}   

y.weekifr = aggregate(weekgroup$x, by=list(week=weekgroup$week), FUN=calcifr.y, drop=F)
o.weekifr = aggregate(weekgroup$x, by=list(week=weekgroup$week), FUN=calcifr.o, drop=F)
all.weekifr = aggregate(weekgroup$x, by=list(week=weekgroup$week), FUN=calcifr, drop=F)
weekifr <- data.frame(y.weekifr$week, y.weekifr$x, o.weekifr$x, all.weekifr$x)
weekifr$index = (seq(3,length(o.dmorti)*2,7))[1:(length(y.weekifr$week))]

m.yifr <- smooth.spline(weekifr$index, weekifr$y.weekifr.x, df=5)
yf <- data.frame(index=seq(1:length(o.dmorti)))
y.ifr <- as.numeric(unlist(predict(m.yifr, yf)$y))

m.oifr <- smooth.spline(weekifr$index, weekifr$o.weekifr.x, df=4)
o.ifr <- as.numeric(unlist(predict(m.oifr, yf)$y))

m.ifr <- smooth.spline(weekifr$index, weekifr$all.weekifr.x, df=4)
all.ifr <- as.numeric(unlist(predict(m.ifr, yf)$y))

pdf("ifr.pdf", width=15, height=5)
par(mfrow=c(1,3))
x <-seq(dstartdate, dstartdate+length(y.ifr)-1, by=1)
plot(x, y.ifr * 100, type='l', main="Time profile of COVID-19 IFR Sweden (<60y)",
        xlab="Date", ylab="Infection Fatality Rate (%)", ylim=c(0, 0.1))
points(dstartdate + weekifr$index, weekifr$y.weekifr.x * 100)

## o.ifr[1:length(slope)] = o.ifr[1:length(slope)] * slope
## o.ifr[(length(slope) + 1):length(o.ifr)] = o.ifr[(length(slope) + 1):length(o.ifr)] * (1 - ifr.reduction)

plot(x, o.ifr * 100, type='l', main="Time profile of COVID-19 IFR Sweden (>60y)",
     xlab="Date", ylab="Infection Fatality Rate (%)", ylim=c(0, 7))
points(dstartdate + weekifr$index, weekifr$o.weekifr.x * 100)

plot(x, all.ifr * 100, type='l', main="Time profile of COVID-19 IFR Sweden (log scale)",
     xlab="Date", ylab="Infection Fatality Rate (%)", ylim=c(0.02, 2), log="y")
points(dstartdate + weekifr$index, weekifr$all.weekifr.x * 100)

print(paste(c("Current IFR (%): ", all.ifr[length(all.ifr)-1] * 100)))

dev.off()

g.N <- aggregate(s$pop_both, by=list(group=s$group), FUN=min, drop=F)$x
g.aggr <- aggregate(s$death, by=list(group=s$group), FUN=sum, drop=F)
g.deaths <- g.aggr$x
print(g.deaths)
g.infected <- g.deaths / g.ifr
g.pct.infected <- g.infected / g.N * 100

y.N <- sum(g.N[1:6])
o.N <- sum(g.N[7:9])

barplot(g.pct.infected ~ g.aggr$group)
print(sum(g.infected[1:6]) / y.N)
print(sum(g.infected[7:9]) / o.N)

print(dstartdate)

## Date of lockdown phase (1)
lockdown_offset <- as.numeric(as.Date("2020/3/11") - dstartdate)

## over how many days the lockdown is estimated to have occurred
lockdown_transition_period <- 7

print(dmort)

## how many deaths at date of lockdown
total_deaths_at_lockdown <- dmort[max(1, lockdown_offset)]

d3 <- as.numeric(as.Date("2020/6/6") - dstartdate)
d4 <- as.numeric(as.Date("2020/7/6") - dstartdate)
