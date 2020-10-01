library(bayesplot)
library(ggplot2)
library(scales)
require(mcmcse)
require(bayestestR)
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

est.ifr.1 <- function(state, params) {
    y.i = -diff(state$y.S)
    y.l <- round(params[14] + rnorm(1, 0, 2))
    o.i = -diff(state$o.S)
    o.l <- round(params[19] + rnorm(1, 0, 2)) 

    (c(state$y.deadi[(y.l+1):length(y.i)], rep(0, y.l)) + c(state$o.deadi[(o.l+1):length(o.i)], rep(0, o.l))) / (y.i + o.i) * 100
}

est.ifr <- function(state, params) {
    fyifr = params[37]
    
    y.i = -diff(state$y.S)
    L = length(y.i)
    gyifr <- y.ifr * (1 + fyifr * g)
    t1 <- state$offset
    t2 <- min(L, t1 + length(gyifr) - 1)
    ydead <- rep(0, L)
    ydead[1:t1] = y.i[1:t1] * gyifr[1]
    ydead[t1:t2] = y.i[t1:t2] * gyifr[1:(t2 - t1 + 1)]
    ydead[t2:L] = y.i[t2:L] * gyifr[length(gyifr)]

    o.i = -diff(state$o.S)
    L = length(o.i)
    goifr <- o.ifr
    t1 <- state$offset
    t2 <- min(L, t1 + length(goifr) - 1)
    odead <- rep(0, L)
    odead[1:t1] = o.i[1:t1] * goifr[1]
    odead[t1:t2] = o.i[t1:t2] * goifr[1:(t2 - t1 + 1)]
    odead[t2:L] = o.i[t2:L] * goifr[length(goifr)]

    (ydead + odead) / (y.i + o.i) * 100
}

est.age <- function(state, params) {
    ifr <- est.ifr(state, params)
    age <- agepop.from.ifr(ifr)
    age
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

all_plots <- function(date_markers) {
    p1 <- makePlot(data_sample, c(dstartdate, plot_end_date),
                   function(state, params) state$deadi, "#3366FF",
                   c("Count", "Deaths"), date_markers, NULL)
    dm1 <- numeric(length(p1$data$x))
    dm1[1:length(p1$data$x)] = NA
    dm1[1:length(dmort)] = dmorti
    p1 <- p1 + geom_line(aes(y = dm1)) + geom_point(aes(y = dm1), size=0.5, color="#1144CC")
    
    p2 <- makePlot(data_sample, c(dstartdate, plot_end_date),
                   function(state, params) state$died, "#3366FF",
                   c("Count", "Total deaths"), date_markers, NULL)
    dm2 <- numeric(length(p2$data$x))
    dm2[1:length(p2$data$x)] = NA
    dm2[1:length(dmort)] = dmort
    p2 <- p2 + geom_line(aes(y = dm2)) + geom_point(aes(y = dm2),  size=0.5, color="#1144CC")

    p3 <- makePlot(data_sample, c(dstartdate, plot_end_date),
                   function(state, params) state$casei, "#FF6633",
                   c("Count", paste(c(CaseLabel, "per day"))), date_markers, NULL)
    dm3 <- numeric(length(p3$data$x))
    dm3[1:length(p3$data$x)] = NA
    dm3[1:length(dcasei)] = dcasei
    p3 <- p3 + geom_line(aes(y = dm3)) + geom_point(aes(y = dm3),  size=0.5, color="#CC4411")

    p4 <- makePlot(data_sample, c(dstartdate, plot_end_date),
                   function(state, params) { (state$E + state$In + state$Is)/N * 100 },
                   "#FFFF66", c(paste(country_adjective, "population (%)"), "Infected people"), date_markers, NULL)

    ## p4 <- makePlot(data_sample, c(dstartdate, plot_end_date),
    ##                function(state, params) { (state$E + state$I)/N * 100 },
    ##                "#FFFF66", c(paste(country_adjective, "population (%)"), "Infected people"), date_markers, NULL)

    p5 <- makePlot(data_sample, c(dstartdate, plot_end_date),
                   function(state, params) { state$R/N * 100 }, "#33FF66",
                   c(paste(country_adjective, "population (%)"), "Recovered from disease"), date_markers, NULL)

    p6 <- makePlot(data_sample, c(dstartdate, plot_end_date),
                   function(state, params) { state$Re }, "#33FFFF",
                   c("Re", "Effective reproduction number (Re)"), date_markers, NULL)

    p6 <- p6 + coord_cartesian(ylim = c(0, NA)) +
        geom_hline(yintercept=1, linetype="solid", color="gray", size=0.5)

    ## p7 <- makePlot(data_sample, c(dstartdate, plot_end_date),
    ##                function(state, params) { state$BS }, "#33FF66",
    ##                c(paste(country_adjective, "BS"), ""), date_markers, NULL)

    ## p8 <- makePlot(data_sample, c(dstartdate, plot_end_date),
    ##                function(state, params) { state$Ib * state$BS }, "#33FF66",
    ##                c(paste(country_adjective, "Ib * BS"), ""), date_markers, NULL)

    ## p9 <- makePlot(data_sample, c(dstartdate, plot_end_date),
    ##                function(state, params) { (state$Ib / (N / state$BS)) }, "#33FF66",
    ##                c("homogeneity", ""), date_markers, NULL)
   
    ## grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow=3)
    grid.arrange(p1, p2, p3, p4, p5, p6, nrow=3)
}

all_plots_age <- function(date_markers) {
    dateRange <- c(plot_start_date, plot_end_date)

    ageGroupLabels = c('< 65y','>= 65y','All (50%, 90% cri)')

    ## Plot 1 : daily deaths

    p1 <- makePlot(data_sample, dateRange,
                   function(state, params) { state$y.deadi + state$o.deadi }, "#3366FF",
                   c("Count", "Deaths"), date_markers, 'solid')

    start <- as.numeric(dstartdate - dateRange[1])

    ##
    ## dm1 <- numeric(length(p1$data$x))
    ## dm1[1:length(p1$data$x)] = NA
    ## dm1[(start + 1):(start + length(dmorti))] = dmorti
    ## p1 <- p1 + geom_line(aes(y = dm1)) + geom_point(aes(y = dm1),  size=0.5, color="#1144CC")

    if (zoom == 1) {
        p1 <- p1 + coord_cartesian(ylim = c(0, 50))
    } else if (zoom == 2) {
        p1 <- p1 + coord_cartesian(xlim = c(as.Date("2020/6/1"), plot_end_date), ylim = c(0, 50))
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
    
    ## ## Plot 2 : total deaths
    
    ## p2 <- makePlot(data_sample, dateRange,
    ##                function(state, params) { state$y.died + state$o.died }, "#3366FF",
    ##                c("Count", "Total deaths"), date_markers, 'solid')
    ## p2 <- p2 + theme(legend.position = c(0.2, 0.85))

    ## dm2 <- numeric(length(p2$data$x))
    ## dm2[1:length(p2$data$x)] = NA
    ## dm2[(start + 1):(start + length(dmort))] = dmort
    ## p2 <- p2 + geom_line(aes(y = dm2)) + geom_point(aes(y = dm2), size=0.25, color="#1144CC")

    ## ## Add y curves
    ## p2 <- addExtraPlot(p2, data_sample, dateRange, function(state, params) { state$y.died }, "#3366FF", "dashed")

    ## dm2y <- numeric(length(p2$data$x))
    ## dm2y[1:length(p2$data$x)] = NA
    ## dm2y[(start + 1):(start + length(y.dmort))] = y.dmort
    
    ## p2 <- p2 + geom_line(aes(y = dm2y), linetype="dashed") + geom_point(aes(y = dm2y),  size=0.25, color="#1144CC")

    ## ## Add o curves
    ## p2 <- addExtraPlot(p2, data_sample, dateRange, function(state, params) { state$o.died }, "#3366FF", "dotted")

    ## dm2o <- numeric(length(p2$data$x))
    ## dm2o[1:length(p2$data$x)] = NA
    ## dm2o[(start + 1):(start + length(o.dmort))] = o.dmort

    ## p2 <- p2 + geom_line(aes(y = dm2o), linetype="dotted") + geom_point(aes(y = dm2o),  size=0.25, color="#1144CC")

    ## p2 <- p2 + scale_colour_manual(values = c("#1144CC"), guide=FALSE) +
    ##     scale_linetype_manual(name = 'Age group',
    ##                           values=c('solid'='solid','dashed'='dashed','dotted'='dotted'),
    ##                           labels = ageGroupLabels)

    ## Plot 3 : cases

    p3 <- makePlot(data_sample, dateRange,
                   function(state, params) state$y.casei + state$o.casei, "#FF6633",
                   c("Count", paste(c("Cases"))), date_markers, 'solid')

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

    if (zoom == 1) {
        p3 <- p3 + coord_cartesian(ylim = c(0, 4000))
    } else if (zoom == 2) {
        p3 <- p3 + coord_cartesian(xlim = c(as.Date("2020/6/1"), plot_end_date),
                                   ylim = c(0, 4000))
    }


    # p3 <- p3 + scale_y_continuous(limits=c(0, 2500))

    ## Plot 3b : daily hospitalisations

    p3b <- makePlot(data_sample, dateRange,
                    function(state, params) { state$y.hospi + state$o.hospi }, "#581845",
                    c("Count", "Hospitalisations"), date_markers, 'solid')

    start <- as.numeric(dstartdate - dateRange[1])
   
    dm3b <- numeric(length(p3b$data$x))
    dm3b[1:length(p3b$data$x)] = NA
    dm3b[(start + 1):(start + length(dhospi))] = dhospi
    p3b <- p3b + geom_line(aes(y = dm3b), size=0.1) + geom_point(aes(y = dm3b),  size=0.25, color="#581845")

    if (zoom == 1) {
        p3b <- p3b + coord_cartesian(ylim = c(0, 300))
    } else if (zoom == 2) {
        p3b <- p3b + coord_cartesian(xlim = c(as.Date("2020/6/1"), plot_end_date), ylim = c(0, 200))
    }

    ## Add y/o curves
    p3b <- addExtraPlot(p3b, data_sample, dateRange, function(state, params) { state$y.hospi }, "#581845", "dashed")
    p3b <- addExtraPlot(p3b, data_sample, dateRange, function(state, params) { state$o.hospi }, "#581845", "dotted")

    p3b <- p3b + scale_colour_manual(values = c("#581845"), guide=FALSE) +
        scale_linetype_manual(name = 'Age group',
                              values=c('solid'='solid','dashed'='dashed','dotted'='dotted'),
                              labels = ageGroupLabels)
    
    ## Plot 4 : infected people
    
    p4 <- makePlot(data_sample, dateRange,
                   function(state, params) { (state$y.E + state$y.I + state$o.E + state$o.I)/N * 100 },
                   "#A67514", c("Population group (%)", "Infected people"), date_markers, 'solid')

    if (zoom == 1) {
        p4 <- p4 + coord_cartesian(ylim = c(0, 3))
    } else if (zoom == 2) {
        p4 <- p4 + coord_cartesian(xlim = c(as.Date("2020/6/1"), plot_end_date),
                                   ylim = c(0, 3))
    }


    ## Add y/o curves
    p4 <- addExtraPlot(p4, data_sample, dateRange, function(state, params) { (state$y.E + state$y.I)/y.N * 100 }, "#A67514", "dashed")
    p4 <- addExtraPlot(p4, data_sample, dateRange, function(state, params) { (state$o.E + state$o.I)/o.N * 100 }, "#A67514", "dotted")

    p4 <- p4 + scale_colour_manual(values = c("#A67514"), guide=FALSE) +
        scale_linetype_manual(name = 'Age group',
                              values=c('solid'='solid','dashed'='dashed','dotted'='dotted'),
                              labels = ageGroupLabels)
    
    ## Plot 5 : recovered people

    p5 <- makePlot(data_sample, dateRange,
                   function(state, params) { (state$y.R + state$o.R)/N * 100 }, "#33CC66",
                   c("Population group (%)", "Recovered from disease"), date_markers, 'solid')
    p5 <- p5 + theme(legend.position = c(0.2, 0.85))

    ## Add y/o curves
    p5 <- addExtraPlot(p5, data_sample, dateRange, function(state, params) { state$y.R/y.N * 100 }, "#33CC66", "dashed")
    p5 <- addExtraPlot(p5, data_sample, dateRange, function(state, params) { state$o.R/o.N * 100 }, "#33CC66", "dotted")

    p5 <- p5 + scale_colour_manual(values = c("#33CC66"), guide=FALSE) +
        scale_linetype_manual(name = 'Age group',
                              values=c('solid'='solid','dashed'='dashed','dotted'='dotted'),
                              labels = ageGroupLabels)

    p6 <- makePlot(data_sample, dateRange,
                   function(state, params) { state$Re }, "#33FFFF",
                   c("Re", "Effective reproduction number (Re)"), date_markers, NULL)

    if (zoom == 1) {
        p6 <- p6 + coord_cartesian(ylim = c(0, 2))
    } else if (zoom == 2) {
        p6 <- p6 + coord_cartesian(xlim = c(as.Date("2020/6/1"), plot_end_date),
                                   ylim = c(0, 2))
    } else {
        p6 <- p6 + coord_cartesian(ylim = c(0, NA))
    }

    p6 <- p6 +  geom_hline(yintercept=1, linetype="solid", color="gray", size=0.5)

    pifr <- makePlot(data_sample, dateRange, est.ifr,
                     "#333333", c("IFR (%)", "Infection fatality rate of infected population"), dates, 'solid')
    pifr <- pifr + theme(legend.position = "none")

    page <- makePlot(data_sample, dateRange, est.age,
                     "#66CC66", c("Mean age", "Mean age of infected population"), dates, 'solid')
    page <- page + theme(legend.position = "none")

    if (zoom == 2) {
        pifr <- pifr + coord_cartesian(xlim = c(as.Date("2020/6/1"), plot_end_date),
                                       ylim = c(0, 0.5))
        page <- page + coord_cartesian(xlim = c(as.Date("2020/6/1"), plot_end_date),
                                       ylim = c(0, 50))
    } else {
        pifr <- pifr + coord_cartesian(ylim = c(0, 1.5))
        page <- page + coord_cartesian(ylim = c(0, 50))
    }
    
    ##grid.arrange(p1, p2, p3, p4, p5, p6, nrow=3)
    grid.arrange(p1, p3b, p3, p4, p6, p5, pifr, page, nrow=2)
}

death_hosp_plots_age <- function(date_markers) {
    dateRange <- c(plot_start_date, plot_end_date)

    ageGroupLabels = c('< 65y','>= 65y','All (50%, 90% cri)')

    ## Plot 1 : daily deaths

    p1 <- makePlot(data_sample, dateRange,
                   function(state, params) { state$y.deadi + state$o.deadi }, "#3366FF",
                   c("Count", "Deaths"), date_markers, 'solid')

    start <- as.numeric(dstartdate - dateRange[1])

    if (zoom == 1) {
        p1 <- p1 + coord_cartesian(ylim = c(0, 50))
    } else if (zoom == 2) {
        p1 <- p1 + coord_cartesian(xlim = c(as.Date("2020/6/1"), plot_end_date), ylim = c(0, 50))
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
    
    ## Plot 3b : daily hospitalisations

    p3b <- makePlot(data_sample, dateRange,
                    function(state, params) { state$y.hospi + state$o.hospi }, "#581845",
                    c("Count", "Hospitalisations"), date_markers, 'solid')

    start <- as.numeric(dstartdate - dateRange[1])
   
    dm3b <- numeric(length(p3b$data$x))
    dm3b[1:length(p3b$data$x)] = NA
    dm3b[(start + 1):(start + length(dhospi))] = dhospi
    p3b <- p3b + geom_line(aes(y = dm3b), size=0.1) + geom_point(aes(y = dm3b),  size=0.25, color="#581845")

    if (zoom == 1) {
        p3b <- p3b + coord_cartesian(ylim = c(0, 300))
    } else if (zoom == 2) {
        p3b <- p3b + coord_cartesian(xlim = c(as.Date("2020/6/1"), plot_end_date), ylim = c(0, 200))
    }

    ## Add y/o curves
    p3b <- addExtraPlot(p3b, data_sample, dateRange, function(state, params) { state$y.hospi }, "#581845", "dashed")
    p3b <- addExtraPlot(p3b, data_sample, dateRange, function(state, params) { state$o.hospi }, "#581845", "dotted")

    p3b <- p3b + scale_colour_manual(values = c("#581845"), guide=FALSE) +
        scale_linetype_manual(name = 'Age group',
                              values=c('solid'='solid','dashed'='dashed','dotted'='dotted'),
                              labels = ageGroupLabels)
    
    ## Plot 4 : infected people
    
    p4 <- makePlot(data_sample, dateRange,
                   function(state, params) { (state$y.E + state$y.I + state$o.E + state$o.I)/N * 100 },
                   "#A67514", c("Population group (%)", "Infected people"), date_markers, 'solid')

    if (zoom == 1) {
        p4 <- p4 + coord_cartesian(ylim = c(0, 2))
    } else if (zoom == 2) {
        p4 <- p4 + coord_cartesian(xlim = c(as.Date("2020/6/1"), plot_end_date),
                                   ylim = c(0, 2))
    }


    ## Add y/o curves
    p4 <- addExtraPlot(p4, data_sample, dateRange, function(state, params) { (state$y.E + state$y.I)/y.N * 100 }, "#A67514", "dashed")
    p4 <- addExtraPlot(p4, data_sample, dateRange, function(state, params) { (state$o.E + state$o.I)/o.N * 100 }, "#A67514", "dotted")

    p4 <- p4 + scale_colour_manual(values = c("#A67514"), guide=FALSE) +
        scale_linetype_manual(name = 'Age group',
                              values=c('solid'='solid','dashed'='dashed','dotted'='dotted'),
                              labels = ageGroupLabels)
    
    ## Plot 5 : recovered people

    p5 <- makePlot(data_sample, dateRange,
                   function(state, params) { (state$y.R + state$o.R)/N * 100 }, "#33CC66",
                   c("Population group (%)", "Recovered from disease"), date_markers, 'solid')
    p5 <- p5 + theme(legend.position = c(0.2, 0.85))

    ## Add y/o curves
    p5 <- addExtraPlot(p5, data_sample, dateRange, function(state, params) { state$y.R/y.N * 100 }, "#33CC66", "dashed")
    p5 <- addExtraPlot(p5, data_sample, dateRange, function(state, params) { state$o.R/o.N * 100 }, "#33CC66", "dotted")

    p5 <- p5 + scale_colour_manual(values = c("#33CC66"), guide=FALSE) +
        scale_linetype_manual(name = 'Age group',
                              values=c('solid'='solid','dashed'='dashed','dotted'='dotted'),
                              labels = ageGroupLabels)

    p6 <- makePlot(data_sample, dateRange,
                   function(state, params) { state$Re }, "#33FFFF",
                   c("Re", "Effective reproduction number (Re)"), date_markers, NULL)

    if (zoom == 1) {
        p6 <- p6 + coord_cartesian(ylim = c(0, 2))
    } else if (zoom == 2) {
        p6 <- p6 + coord_cartesian(xlim = c(as.Date("2020/6/1"), plot_end_date),
                                   ylim = c(0, 2))
    } else {
        p6 <- p6 + coord_cartesian(ylim = c(0, NA))
    }

    p6 <- p6 +  geom_hline(yintercept=1, linetype="solid", color="gray", size=0.5)

    
    grid.arrange(p1, p3b, p4, p6, p5, nrow=2)
}

all_plots_age3 <- function(dates) {
    dateRange <- c(plot_start_date, plot_end_date)

    ageGroupLabels = c('< 40y', '40 - 60y','> 60y','All (50%, 90% cri)')

    ## Plot 1 : daily deaths

    p1 <- makePlot(data_sample, dateRange,
                   function(state, params) { state$y.deadi + state$m.deadi + state$o.deadi }, "#3366FF",
                   c("Count", "Deaths"), dates, 'solid')

    start <- as.numeric(dstartdate - dateRange[1])
   
    dm1 <- numeric(length(p1$data$x))
    dm1[1:length(p1$data$x)] = NA
    dm1[(start + 1):(start + length(dmorti))] = dmorti
    p1 <- p1 + geom_line(aes(y = dm1)) + geom_point(aes(y = dm1),  size=0.5, color="#1144CC")

    if (zoom > 0) {
        p1 <- p1 + coord_cartesian(ylim = c(0, 50))
    }

    ## Add y curves
    p1 <- addExtraPlot(p1, data_sample, dateRange, function(state, params) { state$y.deadi }, "#3366FF", "dashed")

    dm1y <- numeric(length(p1$data$x))
    dm1y[1:length(p1$data$x)] = NA
    dm1y[(start + 1):(start + length(y.dmorti))] = y.dmorti
    
    p1 <- p1 + geom_line(aes(y = dm1y), linetype="dashed") + geom_point(aes(y = dm1y),  size=0.5, color="#1144CC")

    ## Add m curves
    p1 <- addExtraPlot(p1, data_sample, dateRange, function(state, params) { state$m.deadi }, "#3366FF", "dotdash")

    dm1m <- numeric(length(p1$data$x))
    dm1m[1:length(p1$data$x)] = NA
    dm1m[(start + 1):(start + length(m.dmorti))] = m.dmorti
    
    p1 <- p1 + geom_line(aes(y = dm1m), linetype="dotdash") + geom_point(aes(y = dm1m),  size=0.5, color="#1144CC")

    ## Add o curves
    p1 <- addExtraPlot(p1, data_sample, dateRange, function(state, params) { state$o.deadi }, "#3366FF", "dotted")

    dm1o <- numeric(length(p1$data$x))
    dm1o[1:length(p1$data$x)] = NA
    dm1o[(start + 1):(start + length(o.dmorti))] = o.dmorti

    p1 <- p1 + geom_line(aes(y = dm1o), linetype="dotted") + geom_point(aes(y = dm1o),  size=0.5, color="#1144CC")

    p1 <- p1 + scale_colour_manual(values = c("#3366FF"), guide=FALSE) +
        scale_linetype_manual(name = 'Age group',
                              values=c('solid'='solid','dashed'='dashed','dotdash'='dotdash','dotted'='dotted'),
                              labels = ageGroupLabels)
    ## Plot 2 : total deaths
    
    p2 <- makePlot(data_sample, dateRange,
                   function(state, params) { state$y.died + state$m.died + state$o.died }, "#3366FF",
                   c("Count", "Total deaths"), dates, 'solid')
    p2 <- p2 + theme(legend.position = c(0.2, 0.85))

    dm2 <- numeric(length(p2$data$x))
    dm2[1:length(p2$data$x)] = NA
    dm2[(start + 1):(start + length(dmort))] = dmort
    p2 <- p2 + geom_line(aes(y = dm2)) + geom_point(aes(y = dm2), size=0.5, color="#1144CC")

    ## Add y curves
    p2 <- addExtraPlot(p2, data_sample, dateRange, function(state, params) { state$y.died }, "#3366FF", "dashed")

    dm2y <- numeric(length(p2$data$x))
    dm2y[1:length(p2$data$x)] = NA
    dm2y[(start + 1):(start + length(y.dmort))] = y.dmort
    
    p2 <- p2 + geom_line(aes(y = dm2y), linetype="dashed") + geom_point(aes(y = dm2y),  size=0.5, color="#1144CC")

    ## Add m curves
    p2 <- addExtraPlot(p2, data_sample, dateRange, function(state, params) { state$m.died }, "#3366FF", "dotdash")

    dm2 <- numeric(length(p2$data$x))
    dm2[1:length(p2$data$x)] = NA
    dm2[(start + 1):(start + length(m.dmort))] = m.dmort
    
    p2 <- p2 + geom_line(aes(y = dm2), linetype="dotdash") + geom_point(aes(y = dm2),  size=0.5, color="#1144CC")

    ## Add o curves
    p2 <- addExtraPlot(p2, data_sample, dateRange, function(state, params) { state$o.died }, "#3366FF", "dotted")

    dm2o <- numeric(length(p2$data$x))
    dm2o[1:length(p2$data$x)] = NA
    dm2o[(start + 1):(start + length(o.dmort))] = o.dmort

    p2 <- p2 + geom_line(aes(y = dm2o), linetype="dotted") + geom_point(aes(y = dm2o),  size=0.5, color="#1144CC")

    p2 <- p2 + scale_colour_manual(values = c("#1144CC"), guide=FALSE) +
        scale_linetype_manual(name = 'Age group',
                              values=c('solid'='solid','dashed'='dashed','dotdash'='dotdash','dotted'='dotted'),
                              labels = ageGroupLabels)

    ## Plot 3 : cases

    p3 <- makePlot(data_sample, dateRange,
                   function(state, params) { state$y.casei + state$m.casei + state$o.casei }, "#FF6633",
                   c("Count", paste(c("Cases"))), dates, 'solid')
    dm3 <- numeric(length(p3$data$x))
    dm3[1:length(p3$data$x)] = NA
    dm3[(start + 1):(start + length(dcasei))] = dcasei
    p3 <- p3 + geom_line(aes(y = dm3)) + geom_point(aes(y = dm3),  size=0.5, color="#CC4411")

    if (zoom > 0) {
        p3 <- p3 + coord_cartesian(ylim = c(0, 1500))
    }

    ## Add y/m/o curves
    p3 <- addExtraPlot(p3, data_sample, dateRange, function(state, params) { state$y.casei }, "#FF6633", "dashed")
    p3 <- addExtraPlot(p3, data_sample, dateRange, function(state, params) { state$m.casei }, "#FF6633", "dotdash")
    p3 <- addExtraPlot(p3, data_sample, dateRange, function(state, params) { state$o.casei }, "#FF6633", "dotted")

    p3 <- p3 + scale_colour_manual(values = c("#FF6633"), guide=FALSE) +
        scale_linetype_manual(name = 'Age group',
                              values=c('solid'='solid','dashed'='dashed','dotdash'='dotdash','dotted'='dotted'),
                              labels = ageGroupLabels)

    # p3 <- p3 + scale_y_continuous(limits=c(0, 2500))

    ## Plot 4 : infected people
    
    p4 <- makePlot(data_sample, dateRange,
                   function(state, params) { (state$y.E + state$y.I + state$m.E + state$m.I + state$o.E + state$o.I)/N * 100 },
                   "#A67514", c("Population group (%)", "Infected people"), dates, 'solid')

    if (zoom > 0) {
        p4 <- p4 + coord_cartesian(ylim = c(0, 0.5))
    }

    ## Add y/m/o curves
    p4 <- addExtraPlot(p4, data_sample, dateRange, function(state, params) { (state$y.E + state$y.I)/y.N * 100 }, "#A67514", "dashed")
    p4 <- addExtraPlot(p4, data_sample, dateRange, function(state, params) { (state$m.E + state$m.I)/m.N * 100 }, "#A67514", "dotdash")
    p4 <- addExtraPlot(p4, data_sample, dateRange, function(state, params) { (state$o.E + state$o.I)/o.N * 100 }, "#A67514", "dotted")

    p4 <- p4 + scale_colour_manual(values = c("#A67514"), guide=FALSE) +
        scale_linetype_manual(name = 'Age group',
                              values=c('solid'='solid','dashed'='dashed','dotdash'='dotdash','dotted'='dotted'),
                              labels = ageGroupLabels)
    
    ## Plot 5 : recovered people

    p5 <- makePlot(data_sample, dateRange,
                   function(state, params) { (state$y.R + state$m.R + state$o.R)/N * 100 }, "#33CC66",
                   c("Population group (%)", "Recovered from disease"), dates, 'solid')
    p5 <- p5 + theme(legend.position = c(0.2, 0.85))

    ## Add y/m/o curves
    p5 <- addExtraPlot(p5, data_sample, dateRange, function(state, params) { state$y.R/y.N * 100 }, "#33CC66", "dashed")
    p5 <- addExtraPlot(p5, data_sample, dateRange, function(state, params) { state$m.R/m.N * 100 }, "#33CC66", "dotdash")
    p5 <- addExtraPlot(p5, data_sample, dateRange, function(state, params) { state$o.R/o.N * 100 }, "#33CC66", "dotted")

    p5 <- p5 + scale_colour_manual(values = c("#33CC66"), guide=FALSE) +
        scale_linetype_manual(name = 'Age group',
                              values=c('solid'='solid','dashed'='dashed','dotdash'='dotdash','dotted'='dotted'),
                              labels = ageGroupLabels)

    p6 <- makePlot(data_sample, dateRange,
                   function(state, params) { state$Re }, "#33CCCC",
                   c("Re", "Effective reproduction number (Re)"), dates, 'solid')
    ## Add y/m/o curves
    p6 <- addExtraPlot(p6, data_sample, dateRange, function(state, params) { state$y.Re }, "#33CCCC", "dashed")
    p6 <- addExtraPlot(p6, data_sample, dateRange, function(state, params) { state$m.Re }, "#33CCCC", "dotdash")
    p6 <- addExtraPlot(p6, data_sample, dateRange, function(state, params) { state$o.Re }, "#33CCCC", "dotted")
    p6 <- p6 + scale_colour_manual(values = c("#33CCCC"), guide=FALSE) +
        scale_linetype_manual(name = 'Age group',
                              values=c('solid'='solid','dashed'='dashed','dotdash'='dotdash','dotted'='dotted'),
                              labels = ageGroupLabels)

    if (zoom > 0) {
        p6 <- p6 + coord_cartesian(ylim = c(0, 2))
    } else {
        p6 <- p6 + coord_cartesian(ylim = c(0, NA))
    }

    p6 <- p6 +  geom_hline(yintercept=1, linetype="solid", color="gray", size=0.5)

    
    ##grid.arrange(p1, p2, p3, p4, p5, p6, nrow=3)
    grid.arrange(p1, p3, p4, p6, p2, p5, nrow=3)
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

# trimData(4)

quantilePlotSampleSize <- 500
data_sample <- readSample()

## Configure this depending on the model
all_plots <- all_plots_age
##all_plots <- death_hosp_plots_age

## No d6
d6 <- as.numeric(as.Date("2020/12/1") - dstartdate)

plot_start_date <- as.Date("2020/2/13")
plot_end_date <- as.Date("2021/1/1")

date_markers <- data.frame(pos=c(dstartdate + lockdown_offset,
                                 dstartdate + lockdown_offset + lockdown_transition_period,
                                 Sys.Date()),
                           color=c("orange", "red", "darkgray"))

## date_markers <- data.frame(pos=c(dstartdate + lockdown_offset,
##                                  dstartdate + lockdown_offset + lockdown_transition_period,
##                                  dstartdate + d5,
##                                  dstartdate + d7,
##                                  Sys.Date()),
##                            color=c("orange", "red", "orange", "yellow", "darkgray"))
dates <- date_markers

est.Re <- data.frame(quantileData(data_sample, function(state, params) { state$Re }, 0, 250, c(0.05, 0.5, 0.95)))
colnames(est.Re) <- c("q5", "q50", "q95")

print(c(est.Re[Sys.Date() - dstartdate,]))

pdf("current-state-2.pdf", width=25, height=10)

zoom <- 1
all_plots(dates)
zoom <- 0
all_plots(dates)
zoom <- 2
all_plots(dates)

zoom <- 0
plot_end_date <- Sys.Date()

all_plots(dates)

plot_end_date <- as.Date("2021/1/1")

est.died <- data.frame(quantileData(data_sample, function(state, params) { state$y.died + state$o.died }, 0, 400, c(0.05, 0.5, 0.95)))
colnames(est.died) <- c("q5", "q50", "q95")
print(est.died[plot_end_date - dstartdate,])

dmortn <- dmort[length(dmort)]

print(paste("Deaths at date, ", plot_end_date))
print(est.died[plot_end_date - dstartdate,] - dmortn)

## d6s <- c(as.Date("2020/7/29"),
##          as.Date("2020/8/15"),
##          as.Date("2020/9/1"))

## for (i in 1:length(d6s)) {
##     d6d <- d6s[i]
##     d6 <- as.numeric(d6d - dstartdate)

##     dates <- data.frame(pos=c(dstartdate + lockdown_offset,
##                               dstartdate + lockdown_offset + lockdown_transition_period,
##                               Sys.Date(), dstartdate + d6),
##                         color=c("orange", "red", "black", "red"))

##     all_plots(dates)
##     est.died <- data.frame(quantileData(data_sample, function(state, params) { state$y.died + state$o.died }, 0, 300, c(0.05, 0.5, 0.95)))
##     colnames(est.died) <- c("q5", "q50", "q95")
##     print(d6d)
##     print(est.died[plot_end_date - dstartdate,] - dmortn)
## }

dev.off()

## Population group Re values

plot_end_date <- Sys.Date()
dateRange <- c(plot_start_date, plot_end_date)
fair_cols <- c("#38170B","#BF1B0B", "#FFC465", "#66ADE5", "#252A52")

## p7 <- makePlot2(data_sample, dateRange,
##                 function(state, params) { state$y.Re }, fair_cols[2],
##                 c("Re", "Effective reproduction number (Re), < 40y"), dates, 'solid')
## p7 <- p7 + geom_hline(yintercept=1, linetype="solid", color="gray", size=0.5)
## p7 <- addExtraPlotQ(p7, data_sample, dateRange,
##                    function(state, params) { state$m.Re }, fair_cols[3], 'solid')
## p7 <- addExtraPlotQ(p7, data_sample, dateRange,
##                    function(state, params) { state$o.Re }, fair_cols[4], 'solid')
## p7 <- p7 + coord_cartesian(ylim = c(0, 2))
## p7

p7 <- makePlot2(data_sample, dateRange,
                function(state, params) { state$y.Re }, fair_cols[2],
                c("Re", "Effective reproduction number, per age group"), dates, 'dotted')

p7 <- p7 + geom_hline(yintercept=1, linetype="solid", color="gray", size=0.5)
p7 <- addExtraPlotQ(p7, data_sample, dateRange,
                   function(state, params) { state$o.Re }, fair_cols[2], 'dashed')
p7 <- p7 + coord_cartesian(ylim = c(0, 2.5))
p7 <- p7 + scale_linetype_manual(name = "Age group",
                                 values = c("dotted", "dashed"),
                                 labels = c("> 60", "< 60"))
p7 <- p7 + theme(legend.position = c(0.5, 0.9))

pdf("Re-groups.pdf", width=6, height=5)
p7
dev.off()

y.est.infected <- data.frame(quantileData(data_sample, function(state, params) { state$y.E + state$y.I }, 0, lockdown_offset + 250, c(0.05, 0.5, 0.95)))

o.est.infected <- data.frame(quantileData(data_sample, function(state, params) { state$o.E + state$o.I }, 0, lockdown_offset + 250, c(0.05, 0.5, 0.95)))

peaki = 10
print(dstartdate + peaki)
print("Infected at peak 1 young/old:")
print(y.est.infected[10,])
print(o.est.infected[10,])
nowi = as.numeric(Sys.Date() - dstartdate)
print("Infected today young/old:")
print(y.est.infected[nowi,])
print(o.est.infected[nowi,])

est.tr <- function(state, params) {
    ycase_latency <- ocase_latency <- 12
    y.i = -diff(state$y.S)
    o.i = -diff(state$o.S)

    t1 <- state$offset
    L <- length(y.dcasei)

    result <- rep(NA, length(y.i))

    d1 <- (y.dcasei + o.dcasei) /
        (y.i[(t1 - ycase_latency):(t1 + L - ycase_latency - 1)] +
         o.i[(t1 - ocase_latency):(t1 + L - ocase_latency - 1)])

    d <- data.frame(value=d1)
    d$index = seq(1:length(d1))
    md <- smooth.spline(d$index, d$value, df=10)
    ddf <- data.frame(index=seq(1:length(d1)))
    d2 <- as.numeric(unlist(predict(md, ddf)$y))

    result[t1:(t1+L-1)] = d2
    result * 100
}

dateRange <- c(plot_start_date, Sys.Date())

pdf("testing.pdf", width=6, height=4)
ptr <- makePlot(data_sample, dateRange, est.tr,
                "violet", c("Diagnosed infections (%)", "Evolution of diagnosed infections"), dates, 'solid')
ptr <- ptr + theme(legend.position = "none")
ptr <- ptr + coord_cartesian(ylim = c(0, NA))
ptr
dev.off()
