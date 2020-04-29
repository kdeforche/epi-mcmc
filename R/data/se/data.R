####################
##
## Read data into dmort, dmorti, dhosp, dhospi:
##  dmort, dmorti: cumulative and incident deaths per day
##  dhosp, dhospi: cumulative and incident hospitalisations per day
##  dstartdate: first entry
##
####################

se.case <- read.csv("covid19_sweden_new_cases.csv")

print(aggregate(se.case$Total, by=list(date=se.case$Date), FUN=sum))
dhospi <- aggregate(se.case$Total, by=list(date=se.case$Date), FUN=sum)$x

dhospi <- dhospi[27:length(dhospi)]
dstartdate <- as.Date("2020-02-26")

dhosp <- cumsum(dhospi)

se.mort <- read.csv("covid19_sweden_new_deaths.csv")
print(aggregate(se.mort$Total, by=list(date=se.mort$Date), FUN=sum))
dmorti <- aggregate(se.mort$Total, by=list(date=se.mort$Date), FUN=sum)$x

dmorti <- c(rep(0, 14), dmorti)

dmort <- cumsum(dmorti)

death_underreporting_factor <- 1

#####################
## Population size
#####################

N <- 11.2E6

#####################
## Lockdown measures
#####################

## Date of lockdown phase
lockdown_offset <- as.numeric(as.Date("2020/3/10") - dstartdate)

## over how many days the lockdown is estimated to have occurred
## Wikipedia: March 10 (first announcements) to March 24 (restaurant rules) 
lockdown_transition_period <- 14

## how many deaths at date of lockdown
total_deaths_at_lockdown <- 1
