## Based on two files:
##  Global Mobility Report from Google
##  Cases/Deaths Report compiled by ECDC
##   https://www.ecdc.europa.eu/en/publications-data/download-todays-data-geographic-distribution-covid-19-cases-worldwide

## Needs a variable 'country2' which is the two-letter country code (e.g. 'BE')

mob <- read.csv("Global_Mobility_Report.csv")

ecdc <- read.csv("ecdc.csv")

mob_countries <- levels(mob$country_region_code)
ecdc_countries <- levels(ecdc$geoId)

countries <- c(intersect(mob_countries, ecdc_countries), "GB")

good_countries <- c()

for (country2 in countries) {
    c <- subset(mob, country_region_code == country2 & sub_region_1 == '')
    
    s <- c$transit_stations_percent_change_from_baseline + c$workplaces_percent_change_from_baseline
    if (length(which(is.na(s))) > 0) {
        print(c("Skipping", country2, "incomplete google mobility data"))
        next
    }

    plot(s, type='l')

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
    print(dhospi)
    print(dmorti)

    total_deaths <- dmort[length(dmort)]
    print(c(country2, "total deaths:", total_deaths))

    if (total_deaths > 15) {
        good_countries <- c(good_countries, country2)
    }
}

print(paste(good_countries, collapse=' '))
