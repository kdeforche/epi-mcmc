####################
##
## Read data into dmort, dmorti, dhosp, dhospi:
##  dmort, dmorti: cumulative and incident deaths per day
##  dhosp, dhospi: cumulative and incident hospitalisations per day
##  dstartdate: first entry
##
####################

dstartdate <- as.Date("2020/2/20")

es.data <- read.csv("covid-19-ES-CCAA-serieHistorica.csv", stringsAsFactors=F)
es.data$CASOS[is.na(es.data$CASOS)] <- 0
es.data$Fallecidos[is.na(es.data$Fallecidos)] <- 0

es.data$date <- as.Date(es.data$FECHA, "%d/%m/%Y")
print(aggregate(es.data$CASOS, by=list(date=es.data$date), FUN=sum))
dhosp <- aggregate(es.data$CASOS, by=list(date=es.data$date), FUN=sum)$x

dhospi <- numeric(length(dhosp))
for (i in 1:length(dhosp)) {
    dhospi[i] <- if (i==1) 0 else dhosp[i] - dhosp[i-1]
}

print(aggregate(es.data$Fallecidos, by=list(date=es.data$date), FUN=sum))
dmort <- aggregate(es.data$Fallecidos, by=list(date=es.data$date), FUN=sum)$x

dmorti <- numeric(length(dmort))
for (i in 1:length(dmort)) {
    dmorti[i] <- if (i==1) 0 else dmort[i] - dmort[i - 1]
}

#####################
## Population size
#####################

N <- 47E6

#####################
## Lockdown measures
#####################

## Date of lockdown phase (1)
lockdown_offset <- as.numeric(as.Date("2020/3/15") - dstartdate)

## over how many days the lockdown is estimated to have occurred (March 30 staying at home)
lockdown_transition_period <- 15

total_deaths_at_lockdown <- dmort[lockdown_offset]

print(c("total deaths at lockdown", total_deaths_at_lockdown))
