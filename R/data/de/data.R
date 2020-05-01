####################
##
## Read data into dmort, dmorti, dhosp, dhospi:
##  dmort, dmorti: cumulative and incident deaths per day
##  dhosp, dhospi: cumulative and incident hospitalisations per day
##  dstartdate: first entry
##
####################

# Data from wikipedia

dstartdate <- as.Date("2020-03-01")

dhospi <- c(51, 33, 38, 52, 109, 185, 150, 163, 265, 348, 424, 485, 693, 733, 1043, 1174, 1144, 1042, 2801, 2958, 2705, 1948, 4062, 4764, 4118, 4954, 5780, 6294, 3965, 4751, 4615, 5453, 6156, 6174, 6082, 5936, 3677, 3834, 4003, 4974, 5323, 4133, 2821, 2537, 2082, 2486, 2866, 3380, 3609, 2458, 1775, 1785, 2237, 2352, 2337, 2055, 1737, 1018, 1144, 1304, 1478)

dhosp <- cumsum(dhospi)

dmorti <- c(0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 1, 2, 0, 3, 4, 1, 0, 0, 8, 11, 16, 8, 31, 28, 35, 49, 55, 72, 64, 66, 128, 149, 140, 145, 141, 184, 92, 173, 254, 246, 266, 171, 129, 126, 170, 285, 315, 299, 242, 184, 110, 194, 281, 215, 227, 179, 140, 110, 163, 202, 173)

dmort <- cumsum(dmorti)

#####################
## Population size
#####################

N <- 83E6

country_adjective <- "German"

#####################
## Lockdown measures
#####################

## Date of lockdown phase
## Estimated from Google Mobility Report
lockdown_offset <- as.numeric(as.Date("2020/3/13") - dstartdate)

## over how many days the lockdown is estimated to have occurred
lockdown_transition_period <- 10

## how many deaths at date of lockdown
total_deaths_at_lockdown <- 5
