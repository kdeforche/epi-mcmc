## Based on two files:
##  Global Mobility Report from Google
##  Cases/Deaths Report compiled by ECDC
##   https://www.ecdc.europa.eu/en/publications-data/download-todays-data-geographic-distribution-covid-19-cases-worldwide

mob <- read.csv("Global_Mobility_Report.csv")

western <- c('AT','BY','BE','BA','CA','HR','CZ','DK','EE','FI','FR','DE','GR','HU','IE','IL','IT','LV','LT','LU','MD','NL','MK','NO','PL','PT','RO','RS','SK','SI','ES','SE','CH','GB','US')

i = 1
pdf("graphs1.pdf", width=12, height=16)
par(mfrow=c(5, 4))
for (country2 in western) {
    c <- subset(mob, country_region_code == country2 & sub_region_1 == '')
    countryname <- c$country_region[1]
    s <- c$transit_stations_percent_change_from_baseline + c$workplaces_percent_change_from_baseline

    optimpw <- function(par) {
        i1 <- floor(par[1])
        i2 <- floor(par[2])
        a <- par[3]
        b <- par[4]

        if (i1 > i2)
            return(Inf)

        if (i1 < 1 || i1 > length(s) || i2 < 1 || i2 > length(s))
            return(Inf)

        if (country2 == "BY" & i1 < 40)
            return(Inf)

        if (country2 == "IT" & i1 < 10)
            return(Inf)
    
        i1 <- max(2, min(i2, i1))
        i2 <- min(length(s)-1, max(i2, i1))

        e1 <- sum(abs(a - s[1:(i1 - 1)]))
        e2 <- sum(abs(seq(a, b, length.out=(i2 - i1 + 1)) - s[i1:i2]))
        e3 <- sum(abs(b - s[(i2 + 1):length(s)]))

        e1 + e2 + e3
    }

    init <- c(20, 30, 0, -60)
    if (country2 == "BY")
        init <- c(41, 42, 0, -60)

    control <- NULL
    control$maxit <- 100000
    control$parscale <- c(1, 1, 1, 1)
    o <- optim(init, optimpw, control=control)
    o <- optim(o$par, optimpw, control=control)
    print(o$par)

    print(c$date[round(o$par[1:2])])

    plot(as.Date(c$date), s, type='l', xlab="Date", ylab="% transit + % workplaces", main=countryname)
    lxi <- c(1, round(o$par[1]), round(o$par[2]), length(c$date))
    ly <- c(o$par[3], o$par[3], o$par[4], o$par[4])
    lines(as.Date(c$date[lxi]), ly, col='blue', lw=2)

    d1 <- as.Date(c$date[round(o$par[1])])
    d2 <- as.Date(c$date[round(o$par[2])])
    abline(v=d1, col='orange')
    abline(v=d2, col='red')

    mtext(c("D1", "D2"), at=c(d1, d2), side=c(1), cex=0.5)

    y1 <- max(s)
    y2 <- max(s) - 0.1 * (max(s) - min(s))
    
    text(x=max(as.Date(c$date)), y=y1, paste("D1=",substr(d1, 6, 10)), adj=c(1, NA))
    text(x=max(as.Date(c$date)), y=y2, paste("D2=",substr(d2, 6, 10)), adj=c(1, NA))

    i = i + 1

    if (i == 21) {
        dev.off()
        pdf("graphs2.pdf", width=12, height=16)
        par(mfrow=c(5, 4))
    }
}

dev.off()
