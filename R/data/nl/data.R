####################
##
## Read data into dmort, dmorti, dhosp, dhospi:
##  dmort, dmorti: cumulative and incident deaths per day
##  dhosp, dhospi: cumulative and incident hospitalisations per day
##  dstartdate: first entry
## Sources for daily deaths and hospitalizations at
##  https://www.rivm.nl/coronavirus-covid-19/grafieken
## Sources for lockdown measures at
##  https://nl.wikipedia.org/wiki/Coronacrisis_in_Nederland
##
####################

dstartdate <- as.Date("2020/2/27")

nl.hosp <- read.csv("in-ziekenhuis-opgenomen-patienten.csv", sep=";")
dhospi <- nl.hosp$tot.en.met.gisteren
dhosp <- cumsum(dhospi)

nl.mort <- read.csv("overledenen-per-dag.csv", sep=";")
dmorti <- nl.mort$tot.en.met.gisteren
dmort <- cumsum(dmorti)

## remove two latest data point due to delayed reporting
dhospi <- dhospi[1:(length(dhospi) - 2)]
dhosp <- cumsum(dhospi)
dmorti <- dmorti[1:(length(dmorti) - 2)]
dmort <- cumsum(dmorti)

print(dhospi)
print(dhosp)
print (dmorti)
print (dmort)
print (length(dmort))
print (length(dhosp))

print(paste("last day morti: ", dstartdate + length(dmorti) - 1))

#####################
## Population size
#####################

N <- 17.0E6
country_adjective <- "Dutch"

#####################
## Lockdown measures
#####################

## Date of lockdown phase (1)
lockdown_offset <- as.numeric(as.Date("2020/3/10") - dstartdate)

## over how many days the lockdown is estimated to have occurred
lockdown_transition_period <- 5

## how many deaths at date of lockdown
total_deaths_at_lockdown <- dmort[lockdown_offset]
