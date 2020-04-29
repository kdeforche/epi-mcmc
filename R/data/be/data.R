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

print(paste("last day morti: ", dstartdate + length(dmorti) - 1))

death_underreporting_factor <- 1

#####################
## Population size
#####################

N <- 11.5E6
country_adjective <- "Belgian"

#####################
## Lockdown measures
#####################

## Date of lockdown phase (1)
lockdown_offset <- as.numeric(as.Date("2020/3/13") - dstartdate)

## over how many days the lockdown is estimated to have occurred
lockdown_transition_period <- 7

## how many deaths at date of lockdown
total_deaths_at_lockdown <- 5
