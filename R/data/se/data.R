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

dhospi <- c(1, 1, 8, 3, 0, 5, 13, 30, 25, 59, 33, 46, 101, 98, 196, 151, 152, 71, 69, 83, 119, 145, 143, 180, 134, 117, 182, 230, 314, 286, 366, 300, 281, 415, 475, 486, 554, 601, 357, 340, 389, 738, 655, 645, 454, 395, 464, 437, 479, 604, 623, 688, 532, 389, 462, 709, 721, 748, 770, 474, 299, 543, 744) # from wikipedia

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

country_adjective <- "Swedish"

#####################
## Lockdown measures
#####################

## Date of lockdown phase
## Estimated from Google Mobility Report (with bias for Stockholm county)
lockdown_offset <- as.numeric(as.Date("2020/3/09") - dstartdate)

## over how many days the lockdown is estimated to have occurred
## taking into account also Wikipedia restaurant rules (March 24...)
lockdown_transition_period <- 10

## how many deaths at date of lockdown
total_deaths_at_lockdown <- 1
