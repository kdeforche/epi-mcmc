library(bayesplot)
library(ggplot2)
library(scales)
require(mcmcse)
require(bayestestR)
library(grid)
require(gridExtra)

source("settings.R")

source(data, chdir=T)
source(fitmodel, chdir=T)

sourceR <- function(file) {
    source(paste(Rdir, file, sep=""))
}

sourceR("lib/libMCMC.R")
sourceR("lib/libfit.R")

###
## IFR and age profile
###

## Levin et. al. Assessing the Age Specificity of Infection Fatality Rates for COVID-19: Systematic Review, Meta-Analysis, and Public Policy Implications
ifr.from.age <- function(age) {
    exp(-7.53 + 0.119 * age)
}

age.from.ifr <- function(ifr) {
    (log(ifr) + 7.53) / 0.119
}

ifr.from.agepop <- function(age.mean, age.sd) {
    a <- rnorm(1000000, mean=age.mean, sd=age.sd)
    pos.a <- a[a > 0]
    mean(ifr.from.age(a))
}

ifragepops <- function() {
    ifrs <- c()
    ages <- seq(10, 50, 1)
    for (a in ages) {
        ifrs <- c(ifrs, ifr.from.agepop(a, a/2))
    }
    data.frame(age=ages, ifr=ifrs)
}

ifragepop <- ifragepops()
ageifrm <- lm(ifragepop$age ~ log(ifragepop$ifr))

agepop.from.ifr <- function(ifr) {
    a <- ageifrm$coefficients[1]
    b <- ageifrm$coefficients[2]

    a + b * log(ifr)
}

est.ifr <- function(state, params) {
    fyifr = params[37]
    
    y.i = state$y.i
    L = length(y.i)
    gyifr <- y.ifr * (1 + fyifr * g)
    t1 <- state$offset
    t2 <- min(L, t1 + length(gyifr) - 1)
    ydead <- rep(0, L)
    ydead[1:t1] = y.i[1:t1] * gyifr[1]
    if (t2 > t1) {
        ydead[t1:t2] = y.i[t1:t2] * gyifr[1:(t2 - t1 + 1)]
    }
    ydead[t2:L] = y.i[t2:L] * gyifr[length(gyifr)]

    o.i = state$o.i
    L = length(o.i)
    goifr <- o.ifr
    t1 <- state$offset
    t2 <- min(L, t1 + length(goifr) - 1)
    odead <- rep(0, L)
    odead[1:t1] = o.i[1:t1] * goifr[1]
    if (t2 > t1) {
        odead[t1:t2] = o.i[t1:t2] * goifr[1:(t2 - t1 + 1)]
    }
    odead[t2:L] = o.i[t2:L] * goifr[length(goifr)]

    (ydead + odead) / (y.i + o.i) * 100
}

est.age <- function(state, params) {
    ifr <- est.ifr(state, params)
    age <- agepop.from.ifr(ifr)
    age
}

calcLogNormCumProfile <- function(mean, sd)
{
    kbegin = 0
    kend = 60

    result = NULL
    result$kbegin = -kend
    result$kend = -kbegin
    result$values = numeric(result$kend - result$kbegin)
    i = 1
    for (k in kbegin:kend) {
        result$values[i] = plnorm(k-1, mean=log(mean), sd=log(sd))
        i = i + 1
    }

    result$values = 1 - result$values

    result
}

##inhospprof <- calcLogNormCumProfile(8.21, 2.44)
inhospprof <- calcLogNormCumProfile(7.2, 2.3)
delta <- 2

est.beds <- function(state, params) {
    s <- convolute(1.15 * (state$y.hospi + state$o.hospi), state$offset, length(state$y.hospi), inhospprof)
    c(rep(0, state$offset - delta), s, rep(NA, delta))
}

est.icubeds <- function(state, params) {
    beds <- est.beds(state, params)

    ## change to 3 and 3.5 -> -25%
    ## with a logis ?
    s1 <- seq(1:length(beds))
    t0 <- as.numeric(state$offset + (as.Date("2020/10/01") - dstartdate))
    logis1 <- (1 - 0.3 / (1 + exp(-(s1 - t0) * 0.070)))
    t0 <- as.numeric(state$offset + (as.Date("2021/03/01") - dstartdate))
    logis2 <- (1 - 0.25 / (1 + exp(-(s1 - t0) * 0.070)))
    
    icu1 <- c(rep(0, 5),
              beds[1:(length(beds) - 5)] / (3.75 * logis2[6:length(logis2)]))
    icu2 <- beds / (4.25 * logis2) * logis1

    (icu1 + icu2) / 2
}

DIC <- function(posterior1) {
  draw = posterior1

  x = draw[0:-7]
  lik = draw$loglikelihood
  lik.fun = function(x) { calclogl(transformParams(x)) }

  D.bar <- -2*mean(lik)
  if(is.vector(x)) theta.bar = mean(x) else theta.bar <- apply(x,2,mean)
  D.hat <- -2*lik.fun(theta.bar)
  pD <- D.bar - D.hat
  pV <- var(-2*lik)/2
  list(DIC=pD+D.bar,DIC2=pV+D.bar,IC=2*pD+D.bar,pD=pD,pV=pV,Dbar=D.bar,Dhat=D.hat)
}

trimData <- function(count) {
    y.dmorti <<- y.dmorti[1:(length(y.dmorti) - count)]
    o.dmorti <<- o.dmorti[1:(length(o.dmorti) - count)]
    dmorti <<- dmorti[1:(length(dmorti) - count)]

    y.dcasei <<- y.dcasei[1:(length(y.dcasei) - count)]
    o.dcasei <<- o.dcasei[1:(length(o.dcasei) - count)]
    dcasei <<- dcasei[1:(length(dcasei) - count)]

    y.dmort <<- y.dmort[1:(length(y.dmort) - count)]
    o.dmort <<- o.dmort[1:(length(o.dmort) - count)]
    dmort <<- dmort[1:(length(dmort) - count)]

    dcase <<- dcase[1:(length(dcase) - count)]
}

source(fitmodel, chdir=T)

sourceR <- function(file) {
    source(paste(Rdir, file, sep=""))
}

sourceR("lib/libMCMC.R")

## Make this for vacation periods

annotateF <- function(plot) {
    plot
    ##+
    ##    annotate("rect", xmin = d2 + dstartdate, xmax = d3 + dstartdate - 30,
    ##             ymin = 0, ymax = Inf, alpha = .05, fill='red') +
    ##    annotate("rect", xmin = dls2, xmax = dle2,
    ##             ymin = 0, ymax = Inf, alpha = .05, fill='red')
}

all_plots_age <- function(date_markers) {
    dateRange <- c(plot_start_date, plot_end_date)

    if (zoom == 0) {
        date_markers1 <- subset(date_markers, color == "black" | color=="red")
    } else {
        date_markers1 <- date_markers
    }

    ageGroupLabels = c('< 65y','>= 65y','All (50%, 90% cri)')
    start <- as.numeric(dstartdate - dateRange[1])

    ## Plot 1 : daily deaths

    p1 <- makePlot(data_sample, dateRange,
                   function(state, params) { state$y.deadi + state$o.deadi }, "#3366FF",
                   c("Count", "Deaths"), date_markers1, 'solid')


    ##
    ## dm1 <- numeric(length(p1$data$x))
    ## dm1[1:length(p1$data$x)] = NA
    ## dm1[(start + 1):(start + length(dmorti))] = dmorti
    ## p1 <- p1 + geom_line(aes(y = dm1)) + geom_point(aes(y = dm1),  size=0.5, color="#1144CC")

    if (zoom == 1) {
        p1 <- p1 + coord_cartesian(ylim = c(0, 200))
        p1 <- p1 + theme(legend.position = c(0.8, 0.85))
    } else if (zoom == 2) {
        p1 <- p1 + coord_cartesian(xlim = c(zoomStartDate, plot_end_date), ylim = c(0, 75))
        p1 <- p1 + theme(legend.position = c(0.2, 0.85))
    } else if (zoom == 3) {
        p1 <- p1 + coord_cartesian(xlim = c(zoomStartDate, plot_end_date), ylim = c(0, 50))
        p1 <- p1 + theme(legend.position = c(0.8, 0.85))
    }

    ## Add y curves
    p1 <- addExtraPlot(p1, data_sample, dateRange, function(state, params) { state$y.deadi }, "#3366FF", "dashed")

    dm1y <- numeric(length(p1$data$x))
    dm1y[1:length(p1$data$x)] = NA
    dm1y[(start + 1):(start + length(y.dmorti))] = y.dmorti
    
    p1 <- p1 + geom_line(aes(y = dm1y), linetype="dashed", size=0.1) + geom_point(aes(y = dm1y), size=0.25, color="#1144CC")

    ## Add o curves
    p1 <- addExtraPlot(p1, data_sample, dateRange, function(state, params) { state$o.deadi }, "#3366FF", "dotted")

    dm1o <- numeric(length(p1$data$x))
    dm1o[1:length(p1$data$x)] = NA
    dm1o[(start + 1):(start + length(o.dmorti))] = o.dmorti

    p1 <- p1 + geom_line(aes(y = dm1o), linetype="dotted", size=0.1) + geom_point(aes(y = dm1o),  size=0.25, color="#1144CC")

    p1 <- p1 + scale_colour_manual(values = c("#3366FF"), guide=FALSE) +
        scale_linetype_manual(name = 'Age group',
                              values=c('solid'='solid','dashed'='dashed','dotted'='dotted'),
                              labels = ageGroupLabels)

    ## Plot 3 : cases

    p3 <- makePlot(data_sample, dateRange,
                   function(state, params) state$y.casei + state$o.casei, "#FF6633",
                   c("Count", paste(c("Cases"))), date_markers1, 'solid')

    ## dm3 <- numeric(length(p3$data$x))
    ## dm3[1:length(p3$data$x)] = NA
    ## dm3[(start + 1):(start + length(dcasei))] = dcasei
    ## p3 <- p3 + geom_line(aes(y = dm3)) + geom_point(aes(y = dm3),  size=0.25, color="#CC4411")

    ## Add y curves
    p3 <- addExtraPlot(p3, data_sample, dateRange, function(state, params) { state$y.casei }, "#FF6633", "dashed")

    dm3y <- numeric(length(p3$data$x))
    dm3y[1:length(p3$data$x)] = NA
    dm3y[(start + 1):(start + length(y.dcasei))] = y.dcasei
    
    p3 <- p3 + geom_line(aes(y = dm3y), linetype="dashed", size=0.1) + geom_point(aes(y = dm3y),  size=0.25, color="#FF6633")

    ## Add o curves
    p3 <- addExtraPlot(p3, data_sample, dateRange, function(state, params) { state$o.casei }, "#FF6633", "dotted")

    dm3o <- numeric(length(p3$data$x))
    dm3o[1:length(p3$data$x)] = NA
    dm3o[(start + 1):(start + length(o.dcasei))] = o.dcasei

    p3 <- p3 + geom_line(aes(y = dm3o), linetype="dotted", size=0.1) + geom_point(aes(y = dm3o),  size=0.25, color="#FF6633")


    p3 <- p3 + scale_colour_manual(values = c("#FF6633"), guide=FALSE) +
        scale_linetype_manual(name = 'Age group',
                              values=c('solid'='solid','dashed'='dashed','dotted'='dotted'),
                              labels = ageGroupLabels)

    if (zoom == 0) {
        p3 <- p3 + theme(legend.position = c(0.2, 0.85))
    } else if (zoom == 2) {
        p3 <- p3 + coord_cartesian(xlim = c(zoomStartDate, plot_end_date),
                                   ylim = c(0, 7500))
        p3 <- p3 + theme(legend.position = c(0.2, 0.85))
    } else if (zoom == 3) {
        p3 <- p3 + coord_cartesian(xlim = c(zoomStartDate, plot_end_date),
                                   ylim = c(0, 5000))
        p3 <- p3 + theme(legend.position = c(0.8, 0.85))
    }


    # p3 <- p3 + scale_y_continuous(limits=c(0, 2500))

    ## Plot 3b : daily hospitalisations

    p3b <- makePlot(data_sample, dateRange,
                    function(state, params) { state$y.hospi + state$o.hospi }, "#581845",
                    c("Count", "Hospitalisations"), date_markers1, 'solid')

    start <- as.numeric(dstartdate - dateRange[1])
   
    dm3b <- numeric(length(p3b$data$x))
    dm3b[1:length(p3b$data$x)] = NA
    dm3b[(start + 1):(start + length(dhospi))] = dhospi
    p3b <- p3b + geom_line(aes(y = dm3b), size=0.1) + geom_point(aes(y = dm3b),  size=0.25, color="#581845")

    if (zoom == 0) {
        p3b <- p3b + coord_cartesian(ylim = c(0, 900))
        p3b <- p3b + theme(legend.position = c(0.2, 0.85))
    } else if (zoom == 2) {
        p3b <- p3b + coord_cartesian(xlim = c(zoomStartDate, plot_end_date), ylim = c(0, 400))
        p3b <- p3b + theme(legend.position = c(0.2, 0.85))
    } else if (zoom == 3) {
        p3b <- p3b + coord_cartesian(xlim = c(zoomStartDate, plot_end_date), ylim = c(0, 150))
        p3b <- p3b + theme(legend.position = c(0.8, 0.85))
    }

    ## Add y/o curves
    p3b <- addExtraPlot(p3b, data_sample, dateRange, function(state, params) { state$y.hospi }, "#581845", "dashed")
    p3b <- addExtraPlot(p3b, data_sample, dateRange, function(state, params) { state$o.hospi }, "#581845", "dotted")

    p3b <- p3b + scale_colour_manual(values = c("#581845"), guide=FALSE) +
        scale_linetype_manual(name = 'Age group',
                              values=c('solid'='solid','dashed'='dashed','dotted'='dotted'),
                              labels = ageGroupLabels)
    
    ## Plot 4 : infected people

    infectedF <- function(state, params) {
        (state$y.E + state$y.I + state$o.E + state$o.I)/N * 100
    }

    p4 <- makePlotAnnotate(data_sample, dateRange, infectedF,
                           "#A67514", c("Population group (%)", "Infected people (model != reality)"),
                           date_markers, 'solid', annotateF)


    if (zoom == 0) {
        p4 <- p4 + theme(legend.position = c(0.2, 0.75))
    } else if (zoom == 2) {
        p4 <- p4 + coord_cartesian(xlim = c(zoomStartDate, plot_end_date),
                                   ylim = c(0, 3))
        p4 <- p4 + theme(legend.position = c(0.2, 0.75))
    } else if (zoom == 3) {
        p4 <- p4 + coord_cartesian(xlim = c(zoomStartDate, plot_end_date),
                                   ylim = c(0, 1))
        p4 <- p4 + theme(legend.position = c(0.8, 0.75))
    }

    ## p4 <- addExtraPlotQ2(p4, data_sample, dateRange,
    ##                      function(state, params) {
    ##                          d <- (state$mt.E + state$mt.I)/N * 100
    ##                          ifelse(d > 0.001, d, NA)
    ##                      },
    ##                      "#581845")

    ## p4 <- addExtraPlotQ2(p4, data_sample, dateRange,
    ##                      function(state, params) {
    ##                          (state$y.E + state$y.I + state$o.E + state$o.I - state$mt.E - state$mt.I)/N * 100
    ##                      },
    ##                      "#185815")

    ## p4 <- p4 + guides(linetype=guide_legend(keywidth = 3, keyheight = 1),
    ##                   colour=guide_legend(keywidth = 3, keyheight = 1))

    ## Add y/o curves
    p4 <- addExtraPlot(p4, data_sample, dateRange, function(state, params) { (state$y.E + state$y.I)/y.N * 100 }, "#A67514", "dashed")
    p4 <- addExtraPlot(p4, data_sample, dateRange, function(state, params) { (state$o.E + state$o.I)/o.N * 100 }, "#A67514", "dotted")

    ## p4 <- p4 + scale_colour_manual(name = 'Strain',
    ##                                values=c("#A67514"="#A67514", "#581845"="#581845", "#185815"="#185815"),
    ##                                labels=c("All", "B.1.1.7", "Other"),
    ##                                breaks=c("#A67514", "#581845", "#185815"))
    p4 <- p4 + scale_colour_manual(values = c("#A67514"), guide=FALSE)    
    p4 <- p4 + scale_linetype_manual(name = 'Age group',
                              values=c('solid'='solid','dashed'='dashed','dotted'='dotted'),
                              labels = ageGroupLabels)
    
    ## Plot 5 : recovered people

    p5 <- makePlot(data_sample, dateRange,
                   function(state, params) { (state$y.R + state$o.R)/N * 100 }, "#33CC66",
                   c("Population group (%)", "Immune (model != reality)"), date_markers1, 'solid')
    p5 <- p5 + theme(legend.position = c(0.2, 0.85))

    ## Add y/o curves
    p5 <- addExtraPlot(p5, data_sample, dateRange, function(state, params) { state$y.R/y.N * 100 }, "#33CC66", "dashed")
    p5 <- addExtraPlot(p5, data_sample, dateRange, function(state, params) { state$o.R/o.N * 100 }, "#33CC66", "dotted")

    p5 <- p5 + scale_colour_manual(values = c("#33CC66"), guide=FALSE) +
        scale_linetype_manual(name = 'Age group',
                              values=c('solid'='solid','dashed'='dashed','dotted'='dotted'),
                              labels = ageGroupLabels)

    if (zoom == 2) {
        p5 <- p5 + coord_cartesian(xlim = c(zoomStartDate, plot_end_date))
    }

    ## Plot 6 : Re
    
    p6 <- makePlotAnnotate(data_sample, dateRange,
                           function(state, params) { state$Re }, "#33FFFF",
                           c("Re", "Effective reproduction number (Re)"), date_markers, NULL,
                           annotateF)

    p6 <- p6 + scale_colour_manual(values = c("#33FFFF"), guide=FALSE)

    ## p6 <- addExtraPlotQ2(p6, data_sample, dateRange,
    ##                      function(state, params) { ifelse(state$mt.Re == 0, NA, state$mt.Re) }, "#581845")

    ## p6 <- addExtraPlotQ2(p6, data_sample, dateRange,
    ##                      function(state, params) { state$wt.Re }, "#185815")

    ## p6 <- p6 + scale_colour_manual(name = 'Strain',
    ##                                values=c("#185815"="#185815", "#33FFFF"="33FFFF", "#581845"="#581845"),
    ##                                labels=c("All", "B.1.1.7", "Other"),
    ##                                breaks=c("#33FFFF", "#581845", "#185815"))
    
    if (zoom == 1) {
        p6 <- p6 + coord_cartesian(ylim = c(0, 2))
    } else if (zoom == 2) {
        p6 <- p6 + coord_cartesian(xlim = c(zoomStartDate, plot_end_date),
                                   ylim = c(0, 1.5))
    } else if (zoom == 3) {
        p6 <- p6 + coord_cartesian(xlim = c(zoomStartDate, plot_end_date),
                                   ylim = c(0, NA))
    } else {
        p6 <- p6 + coord_cartesian(ylim = c(0, NA))
    }

    p6 <- p6 +  geom_hline(yintercept=1, linetype="solid", color="gray", size=0.5)

    p7 <- makePlotAnnotate(data_sample, dateRange,
                           function(state, params) {
                               d <- (state$mt.y.i + state$mt.o.i)/(state$y.i + state$o.i) * 100
                               ifelse(d > 0.001, d, NA)
                           },
                           "#33FFFF",
                           c("Percentage (%)", "Percentage infections by B.1.1.7 variant"), date_markers, NULL,
                           annotateF)
    p7 <- p7 + theme(legend.position = "none")

    if (zoom == 1) {
    } else if (zoom == 2) {
        p7 <- p7 + coord_cartesian(xlim = c(zoomStartDate, plot_end_date),
                                   ylim = c(0, NA))
    }
    
    pifr <- makePlot(data_sample, dateRange, est.ifr,
                     "#333333", c("IFR (%)", "Infection fatality rate of infected population"), dates, 'solid')
    pifr <- pifr + theme(legend.position = "none")

    page <- makePlot(data_sample, dateRange, est.age,
                     "#66CC66", c("Mean age", "Mean age of infected population"), dates, 'solid')
    page <- page + theme(legend.position = "none")


    fair_cols <- c("#38170B","#BF1B0B", "#FFC465", "#237DB5", "#752A22")

    colbeds <- fair_cols[4]
    colicu <- fair_cols[5]

    lbeds <- "Hospital beds"
    licu <- "ICU beds"
    
    pbeds <- makePlot(data_sample, dateRange, est.beds, colbeds,
                      c("Count", "Hospital/ICU bed occupation"), dates, NULL)

    pbeds <- addExtraPlotQ(pbeds, data_sample, dateRange, est.icubeds, colicu, NULL)

    pdbeds <- numeric(length(pbeds$data$x))
    pdbeds[1:length(pbeds$data$x)] = NA
    pdbeds[(start + 1):(start + length(dbeds))] = dbeds
    pbeds <- pbeds + geom_line(aes(y = pdbeds, colour=colbeds), size=0.1) +
        geom_point(aes(y = pdbeds, colour=colbeds), size=0.5)

    pdicu <- numeric(length(pbeds$data$x))
    pdicu[1:length(pbeds$data$x)] = NA
    pdicu[(start + 1):(start + length(dicu))] = dicu
    pbeds <- pbeds + geom_line(aes(y = pdicu, colour=colicu), size=0.1) +
        geom_point(aes(y = pdicu, colour=colicu), size=0.5)

    pbeds <- pbeds + scale_colour_identity(guide = "legend",
                                           labels = c(lbeds, licu)) +
        theme(legend.position = c(0.2, 0.85)) +
        theme(legend.title = element_blank())

    if (zoom == 0) {
        labelDate <- as.Date("2020/2/25")
    } else {
        labelDate <- as.Date("2020/10/7")
    }
    
    ## pbeds <- pbeds + geom_hline(yintercept=1500, linetype="dashed", color=colbeds, size=0.5) +
    ##    geom_text(aes(as.Date("2020/09/01"), 1500, label = "Phase 0", vjust = -1))
    ## pbeds <- pbeds + geom_hline(yintercept=2500, linetype="dashed", color=colbeds, size=0.5) +
    ##     geom_text(aes(as.Date("2020/09/01"), 2500, label = "Phase 1A", vjust = -1))
    pbeds <- pbeds + geom_hline(yintercept=5000, linetype="dashed", color=colbeds, size=0.5) +
        annotate("text", x=labelDate, y=5000, label = "Phase 1B",
                 vjust = -0.3, color=colbeds)
    pbeds <- pbeds + geom_hline(yintercept=7*1500, linetype="dashed", color=colbeds, size=0.5) +
        annotate("text", x=labelDate, y=7*1500, label = "Phase 2A",
                 vjust = -0.3, color=colbeds)
    pbeds <- pbeds + geom_hline(yintercept=7*2000, linetype="dashed", color=colbeds, size=0.5) +
        annotate("text", x=labelDate, y=7*2000, label = "Phase 2B",
                 vjust = -0.3, color=colbeds)

    ## pbeds <- pbeds + geom_hline(yintercept=300, linetype="dashed", color=colicu, size=0.5) +
    ##     geom_text(aes(as.Date("2020/09/01"), 300, label = "Phase 0", vjust = -1))
    ## pbeds <- pbeds + geom_hline(yintercept=500, linetype="dashed", color=colicu, size=0.5) +
    ##     geom_text(aes(as.Date("2020/09/01"), 500, label = "Phase 1A", vjust = -1))
    pbeds <- pbeds + geom_hline(yintercept=1000, linetype="dashed", color=colicu, size=0.5) +
        annotate("text", x=labelDate, y=1000, label = "Phase 1B",
                 vjust = -0.3, color=colicu)
    pbeds <- pbeds + geom_hline(yintercept=1500, linetype="dashed", color=colicu, size=0.5) +
        annotate("text", x=labelDate, y=1500, label = "Phase 2A",
                 vjust = -0.3, color=colicu)
    pbeds <- pbeds + geom_hline(yintercept=2000, linetype="dashed", color=colicu, size=0.5) +
        annotate("text", x=labelDate, y=2000, label = "Phase 2B",
                 vjust = -0.3, color=colicu)

    if (zoom == 0) {
        pbeds <- pbeds + coord_cartesian(ylim = c(0, 10000))
    } else if (zoom == 2) {
        pifr <- pifr + coord_cartesian(xlim = c(zoomStartDate, plot_end_date),
                                       ylim = c(0, 0.5))
        page <- page + coord_cartesian(xlim = c(zoomStartDate, plot_end_date),
                                       ylim = c(0, 50))
        pbeds <- pbeds + coord_cartesian(xlim = c(zoomStartDate, plot_end_date),
                                         ylim = c(0, 4000))
    } else if (zoom == 3) {
        pifr <- pifr + coord_cartesian(xlim = c(zoomStartDate, plot_end_date),
                                       ylim = c(0, 1.0))
        page <- page + coord_cartesian(xlim = c(zoomStartDate, plot_end_date),
                                       ylim = c(0, 50))
        pbeds <- pbeds + coord_cartesian(xlim = c(zoomStartDate, plot_end_date),
                                         ylim = c(0, 1500))
    } else {
        pifr <- pifr + coord_cartesian(ylim = c(0, 1.2))
        page <- page + coord_cartesian(ylim = c(0, 50))
    }

    grid.arrange(p1, p3b, p3, pbeds, p4, p6, p5, pifr, ncol=4,
                 top=textGrob(title, gp=gpar(fontsize=20)))
}

other_plots_age <- function(date_markers) {
    dateRange <- c(plot_start_date, plot_end_date)

    if (zoom == 0) {
        date_markers1 <- subset(date_markers, color == "black" | color=="red")
    } else {
        date_markers1 <- date_markers
    }

    ageGroupLabels = c('< 65y','>= 65y','All (50%, 90% cri)')

    start <- as.numeric(dstartdate - dateRange[1])
    
    ## ## Plot 2 : total deaths
    
    p2 <- makePlot(data_sample, dateRange,
                   function(state, params) { state$y.died + state$o.died }, "#3366FF",
                   c("Count", "Total deaths"), date_markers1, 'solid')
    p2 <- p2 + theme(legend.position = c(0.2, 0.85))

    dm2 <- numeric(length(p2$data$x))
    dm2[1:length(p2$data$x)] = NA
    dm2[(start + 1):(start + length(dmort))] = dmort
    p2 <- p2 + geom_line(aes(y = dm2)) + geom_point(aes(y = dm2), size=0.25, color="#1144CC")

    ## Add y curves
    p2 <- addExtraPlot(p2, data_sample, dateRange, function(state, params) { state$y.died }, "#3366FF", "dashed")

    dm2y <- numeric(length(p2$data$x))
    dm2y[1:length(p2$data$x)] = NA
    dm2y[(start + 1):(start + length(y.dmort))] = y.dmort
    
    p2 <- p2 + geom_line(aes(y = dm2y), linetype="dashed") + geom_point(aes(y = dm2y),  size=0.25, color="#1144CC")

    ## Add o curves
    p2 <- addExtraPlot(p2, data_sample, dateRange, function(state, params) { state$o.died }, "#3366FF", "dotted")

    dm2o <- numeric(length(p2$data$x))
    dm2o[1:length(p2$data$x)] = NA
    dm2o[(start + 1):(start + length(o.dmort))] = o.dmort

    p2 <- p2 + geom_line(aes(y = dm2o), linetype="dotted") + geom_point(aes(y = dm2o),  size=0.25, color="#1144CC")

    p2 <- p2 + scale_colour_manual(values = c("#1144CC"), guide=FALSE) +
        scale_linetype_manual(name = 'Age group',
                              values=c('solid'='solid','dashed'='dashed','dotted'='dotted'),
                              labels = ageGroupLabels)

    p7 <- makePlotAnnotate(data_sample, dateRange,
                           function(state, params) { state$Rt }, "#33FFFF",
                           c("Rt", "Time-varying reproduction number (Rt)"), date_markers, NULL,
                           annotateF)

    p7 <- p7 + scale_colour_manual(values = c("#33FFFF"), guide=FALSE)

    if (zoom == 1) {
        p7 <- p7 + coord_cartesian(ylim = c(0, 2))
    } else if (zoom == 2) {
        p7 <- p7 + coord_cartesian(xlim = c(zoomStartDate, plot_end_date),
                                   ylim = c(0, 2))
    } else {
        p7 <- p7 + coord_cartesian(ylim = c(0, NA))
    }

    p7 <- p7 +  geom_hline(yintercept=1, linetype="solid", color="gray", size=0.5)

    ## Beta

    cols <- c("#00AFBB", "#E7B800", "#FC4E07")

    col.y <- cols[1]
    col.o <- cols[2]
    col.yo <- cols[3]

    pbeta <- makePlot2(data_sample, dateRange, function(state, params) { state$y.beta },
                       col.y, c("Infections per day", "Infections per day per infected individual"),
                       dates, NULL)

    pbeta <- addExtraPlotQ2(pbeta, data_sample, dateRange,
                            function(state, params) { state$o.beta },
                            col.o)

    pbeta <- addExtraPlotQ2(pbeta, data_sample, dateRange,
                            function(state, params) { state$yo.beta },
                            col.yo)

    pbeta <- pbeta + scale_colour_identity(guide = "legend",
                                  labels = c("<65 to <65", ">65 to >65", "<65 to >65")) +
##        theme(legend.position = c(0.2, 0.85)) +
        theme(legend.title = element_blank())
   
    grid.arrange(p2, p7, pbeta, nrow=1)
}

options(scipen=999)

readSample <- function() {
    ## input files: try _1, _2, etc...
    files <- c()
    for (f1 in inputfiles) {
        if (file.exists(f1))
            files <- c(files, f1)
        else {
            for (i in 1:10) {
                f <- paste(f1, "_", i, sep="")
                if (file.exists(f))
                    files <- c(files, f)
            }
        }
    }

    print(files)
   
    all <- readData(files)
    dic = DIC(all)
    print(dic)

    posterior <- all[0:-7]

    plot(ts(subset(posterior, select=keyparamnames)))

    ## compute credibility intervals
    print(data.frame(ci(posterior, ci=0.01)))
    print(data.frame(ci(posterior, ci=0.50)))
    print(data.frame(ci(posterior, ci=0.95)))

    posterior

    selection <- 1:dim(posterior)[1]
    scount <- length(selection)
    draws <- sample(1:scount, min(scount, quantilePlotSampleSize))
    data_sample <- posterior[selection[draws],]

    data_sample
}

outputdir <- "output"

source("settings.R")

quantilePlotSampleSize <- 500
data_sample <- readSample()

all_plots <- all_plots_age

plot_start_date <- as.Date("2020/2/13")

date_markers <- data.frame(pos=c(dstartdate + lockdown_offset,
                                 dstartdate + lockdown_offset + lockdown_transition_period,
                                 dstartdate + d5,
                                 dstartdate + d9,
                                 as.Date("2020/11/1"),
                                 Sys.Date()),
                           color=c("orange", "red", "orange", "orange", "red", "black"))
dates <- date_markers

Range <- 600

est.Re <- data.frame(quantileData(data_sample, function(state, params) { state$Re }, 0, Range, c(0.05, 0.5, 0.95)))

colnames(est.Re) <- c("q5", "q50", "q95")
wt.est.Re <- data.frame(quantileData(data_sample, function(state, params) { state$wt.Re }, 0, Range, c(0.05, 0.5, 0.95)))

mt.est.Re <- data.frame(quantileData(data_sample, function(state, params) { state$mt.Re }, 0, Range, c(0.05, 0.5, 0.95)))

maxRe <- marginalizeData(data_sample, function(state, params) { state$Re }, 0, Range,
                         function(d) {
                             max(d[(Sys.Date() - dstartdate):(Sys.Date() - dstartdate + 60)])
                         })
print("Maximum Re in next 2 months")
print(quantile(maxRe, c(0.05, 0.5, 0.95)))

print("Re Today")
print(unlist(est.Re[Sys.Date() - dstartdate + 1,]))
print("Re @ d10")
print(unlist(est.Re[d10,]))


y.est.infected <- data.frame(quantileData(data_sample, function(state, params) { state$y.E + state$y.I }, 0, lockdown_offset + Range, c(0.05, 0.5, 0.95)))

o.est.infected <- data.frame(quantileData(data_sample, function(state, params) { state$o.E + state$o.I }, 0, lockdown_offset + Range, c(0.05, 0.5, 0.95)))

mt.est.infected <- data.frame(quantileData(data_sample, function(state, params) { state$mt.E + state$mt.I }, 0, lockdown_offset + Range, c(0.05, 0.5, 0.95)))

est.case <- data.frame(quantileData(data_sample, function(state, params) { state$y.casei + state$o.casei }, 0, lockdown_offset + Range, c(0.05, 0.5, 0.95)))

est.hosp <- data.frame(quantileData(data_sample, function(state, params) { state$y.hospi + state$o.hospi }, 0, lockdown_offset + Range, c(0.05, 0.5, 0.95)))

sapply(c(Sys.Date(), as.Date("2021/5/1"), as.Date("2021/6/1"), as.Date("2021/7/1")),
       function(d) {
           do = as.numeric(d - dstartdate)
           ds = as.character(d)
           print(paste(" === Infected young/old/B.1.1.7 on :", ds))
           print(y.est.infected[do,])
           print(o.est.infected[do,])
           print(mt.est.infected[do,])
           print("Re all/wt/B.1.1.7:")
           print(est.Re[do,])
           print(wt.est.Re[do,])
           print(mt.est.Re[do,])
           print("Hosp :")
           print(est.hosp[do,])
           print("Cases :")
           print(est.case[do,])
       })

## est.hosp$ds <- dstartdate + (1:602 - 1)
## ss <- subset(est.hosp, ds > as.Date("2021-03-01"))
## plot(ss$ds, ss$X2)

est.died <- data.frame(quantileData(data_sample, function(state, params) { state$y.died + state$o.died }, 0, Range, c(0.05, 0.5, 0.95)))
colnames(est.died) <- c("q5", "q50", "q95")
print(est.died[plot_end_date - dstartdate,])

dmortn <- dmort[length(dmort)]

print(paste("Deaths at date, ", plot_end_date))
print(est.died[plot_end_date - dstartdate,] - dmortn)

est.totalhosp <- data.frame(quantileData(data_sample, function(state, params) { cumsum(state$y.hosp + state$o.hosp) }, 0, lockdown_offset + Range, c(0.05, 0.5, 0.95)))

print(paste("Total hospitalized at date, ", plot_end_date))
print(est.totalhosp[plot_end_date - dstartdate,])

zoom <- 0

pdf("current-state-2.pdf", width=25, height=10)
##pdf("current-state-2.pdf", width=15, height=10)
all_plots(dates)
zoomStartDate <- as.Date("2021/2/1")
zoom <- 2
##zoom <- 3
all_plots(dates)

dev.off()

pdf("other-plots.pdf", width=19, height=6)
zoom <- 0
other_plots_age(dates)
dev.off()

## Population group Re values

dateRange <- c(plot_start_date, plot_end_date)
fair_cols <- c("#38170B","#BF1B0B", "#FFC465", "#66ADE5", "#252A52")

p7 <- makePlot2(data_sample, dateRange,
                function(state, params) { state$y.Re }, fair_cols[2],
                c("Re", "Effective reproduction number, per age group"), dates, 'dotted')

p7 <- p7 + geom_hline(yintercept=1, linetype="solid", color="gray", size=0.5)
p7 <- addExtraPlotQ(p7, data_sample, dateRange,
                   function(state, params) { state$o.Re }, fair_cols[2], 'dashed')
p7 <- p7 + coord_cartesian(ylim = c(0, 2.5))
p7 <- p7 + scale_linetype_manual(name = "Age group",
                                 values = c("dotted", "dashed"),
                                 labels = c("> 65", "< 65"))
p7 <- p7 + theme(legend.position = c(0.5, 0.9))

pdf("Re-groups.pdf", width=6, height=5)
p7
dev.off()

est.tr <- function(state, params) {
    y.i = state$y.i
    o.i = state$o.i

    ycase_latency <- params[44]
    ocase_latency <- params[45]
    
    ycase_cv_profile = caseProfile(ycase_latency, CLsd)
    ocase_cv_profile = caseProfile(ocase_latency, CLsd)
    padding = max(-ycase_cv_profile$kbegin, -ocase_cv_profile$kbegin) + 1

    y.d = y.i
    y.d[(padding + 1):length(y.i)] = convolute(y.i, padding + 1, length(y.i), ycase_cv_profile)

    o.d = o.i
    o.d[(padding + 1):length(o.i)] = convolute(o.i, padding + 1, length(o.i), ocase_cv_profile)
    
    t1 <- state$offset
    L <- length(y.dcasei)

    result <- rep(NA, length(y.d))

    d1 <- (y.dcasei + o.dcasei) / (y.d[t1:(t1 + L - 1)] + o.d[t1:(t1 + L - 1)])

    d <- data.frame(value=d1)
    d$index = seq(1:length(d1))
    md <- smooth.spline(d$index, d$value, df=15)
    ddf <- data.frame(index=seq(1:length(d1)))
    d2 <- as.numeric(unlist(predict(md, ddf)$y))

    result[t1:(t1+L-1)] = d2
    result * 100
}

dateRange <- c(plot_start_date, Sys.Date())

pdf("testing.pdf", width=6, height=4)
ptr <- makePlot(data_sample, dateRange, est.tr,
                "violet", c("Diagnosed infections (%)", "Evolution of proportion of infections diagnosed using PCR"), dates, 'solid')
ptr <- ptr + theme(legend.position = "none")
ptr <- ptr + coord_cartesian(ylim = c(0, NA))
ptr
dev.off()
