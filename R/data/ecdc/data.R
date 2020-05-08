## Based on two files:
##  Global Mobility Report from Google
##  Cases/Deaths Report compiled by ECDC
##   https://www.ecdc.europa.eu/en/publications-data/download-todays-data-geographic-distribution-covid-19-cases-worldwide

## Needs a variable 'country2' which is the two-letter country code (e.g. 'BE')

mob <- read.csv("Global_Mobility_Report.csv")

c <- subset(mob, country_region_code == country2 & sub_region_1 == '')
s <- c$transit_stations_percent_change_from_baseline + c$workplaces_percent_change_from_baseline
plot(s, type='l')

optimpw <- function(par) {
    i1 <- floor(par[1])
    i2 <- floor(par[2])
    a <- par[3]
    b <- par[4]

    if (i1 > i2)
        return(Inf)
    
    i1 <- max(2, min(i2, i1))
    i2 <- min(length(s)-1, max(i2, i1))

    e1 <- sum(abs(a - s[1:(i1 - 1)]))
    e2 <- sum(abs(seq(a, b, length.out=(i2 - i1 + 1)) - s[i1:i2]))
    e3 <- sum(abs(b - s[(i2 + 1):length(s)]))

    e1 + e2 + e3
}

init <- c(20, 30, 0, -60)

control <- NULL
control$maxit <- 10000
control$parscale <- c(1, 1, 1, 1)
o <- optim(init, optimpw, control=control)
print(o$par)

print(c$date[round(o$par[1:2])])

ecdc <- read.csv("ecdc.csv")

mob_countries <- levels(mob$country_region_code)
ecdc_countries <- levels(ecdc$geoId)

print(intersect(mob_countries, ecdc_countries))

if (country2 == 'GB') {
    country2 = 'UK'
}

country.data <- subset(ecdc, geoId == country2)

country.data$date <- as.Date(paste(country.data$year,country.data$month,country.data$day,sep='/'))

nzcases <- which(country.data$cases > 0)
starti <- nzcases[length(nzcases)]

dstartdate <- country.data$date[starti] ## they are in reversed order

dhospi <- rev(country.data$cases[1:starti])
dmorti <- rev(country.data$deaths[1:starti])

dhosp <- cumsum(dhospi)
dmort <- cumsum(dmorti)

N <- country.data$popData2018[1]
country_adjective <- country.data$countriesAndTerritories[1]

print(dstartdate)

lockdown_offset <- as.numeric(as.Date(c$date[round(o$par[1])]) - dstartdate)
lockdown_transition_period <- round(o$par[2]) - round(o$par[1])
total_deaths_at_lockdown <- max(1, dmort[lockdown_offset])

print(dhospi)
print(dmorti)

print(c("Total deaths:", dmort[length(dmort)]))

if (dmort[length(dmort)] < 15) {
    print("Too few deaths, quiting")
    quit()
}

print(c(lockdown_offset, lockdown_transition_period, total_deaths_at_lockdown))

FitTotalPeriod <- length(dmort) + 90
