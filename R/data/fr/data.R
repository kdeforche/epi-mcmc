####################
##
## Read data into dmort, dmorti, dhosp, dhospi:
##  dmort, dmorti: cumulative and incident deaths per day
##  dhosp, dhospi: cumulative and incident hospitalisations per day
##  dstartdate: first entry
##
## Source : https://www.data.gouv.fr/fr/datasets/donnees-hospitalieres-relatives-a-lepidemie-de-covid-19
##  donnees-hospitalieres-nouveaux-covid19-DATE.csv
## Source population : https://en.wikipedia.org/w/index.php?title=Demographics_of_France&oldid=950365663 
## Lockdown start date : https://en.wikipedia.org/w/index.php?title=2020_coronavirus_pandemic_in_France&oldid=950680978
##
####################

dstartdate <- as.Date("2020/3/19")

fr.hosp <- read.csv("donnees.csv", sep=";")
print(aggregate(fr.hosp$incid_hosp, by=list(date=fr.hosp$jour), FUN=sum))
dhospi <- aggregate(fr.hosp$incid_hosp, by=list(date=fr.hosp$jour), FUN=sum)$x

dhosp <- numeric(length(dhospi))
for (i in 1:length(dhospi)) {
    dhosp[i] <- if (i==1) dhospi[i] else dhosp[i - 1] + dhospi[i]
}

fr.mort <- read.csv("donnees.csv", sep=";")
print(aggregate(fr.mort$incid_dc, by=list(date=fr.mort$jour), FUN=sum))
dmorti <- aggregate(fr.mort$incid_dc, by=list(date=fr.mort$jour), FUN=sum)$x

dmort <- numeric(length(dmorti))
for (i in 1:length(dmort)) {
    dmort[i] <- if (i==1) dmorti[i] else dmort[i - 1] + dmorti[i]
}

#####################
## Population size
#####################

N <- 67.1E6

#####################
## Lockdown measures
#####################

## Date of lockdown phase (1)
lockdown_offset <- as.numeric(as.Date("2020/3/19") - dstartdate)

## over how many days the lockdown is estimated to have occurred
lockdown_transition_period <- 3

total_deaths_at_lockdown <- 372

print (dhospi)
print (dhosp)
print (dmorti)
print (dmort)

country_adjective <- "French"

## how many deaths at date of lockdown
total_deaths_at_lockdown <- 244
