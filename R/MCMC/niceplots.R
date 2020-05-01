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

##
## also plot median line (50% quantile too but very light) for paths not taken
## do not plot age-specific data
## mark lockdown + may 3
## first until July
## then until December
##

bg_color <- "#777777"

all_plots_bg <- function(date_markers) {
    dateRange <- c(dstartdate, plot_end_date)

    p1 <- makePlot(data_sample, c(dstartdate, plot_end_date),
                   function(state) { state$y.deadi + state$o.deadi }, bg_color,
                   c("Count", "Deaths per day"), date_markers, 'dashed')

    p2 <- makePlot(data_sample, c(dstartdate, plot_end_date),
                   function(state) { state$y.died + state$o.died }, bg_color,
                   c("Count", "Total deaths"), date_markers, 'dashed')

    p3 <- makePlot(data_sample, c(dstartdate, plot_end_date),
                   function(state) { state$y.hospi + state$o.hospi }, bg_color,
                   c("Count", paste(c(HospLabel, "per day"))), date_markers, 'dashed')

    p4 <- makePlot(data_sample, c(dstartdate, plot_end_date),
                   function(state) { (state$y.E + state$y.I + state$o.E + state$o.I)/N * 100 },
                   bg_color, c("Belgian population (%)", "Infected people"), date_markers, 'dashed')

    p5 <- makePlot(data_sample, c(dstartdate, plot_end_date),
                   function(state) { (state$y.R + state$o.R)/N * 100 }, bg_color,
                   c("Belgian population (%)", "Recovered from disease"), date_markers, 'dashed')

    list(p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5)
}

add_all_plots_bg <- function(plots) {
    dateRange <- c(dstartdate, plot_end_date)

    p1 <- addExtraPlotQ(plots$p1, data_sample, c(dstartdate, plot_end_date),
                        function(state) { state$y.deadi + state$o.deadi }, bg_color, "solid")

    p2 <- addExtraPlotQ(plots$p2, data_sample, c(dstartdate, plot_end_date),
                        function(state) { state$y.died + state$o.died }, bg_color, "solid")

    p3 <- addExtraPlotQ(plots$p3, data_sample, c(dstartdate, plot_end_date),
                        function(state) { state$y.hospi + state$o.hospi }, bg_color, "solid")

    p4 <- addExtraPlotQ(plots$p4, data_sample, c(dstartdate, plot_end_date),
                        function(state) { (state$y.E + state$y.I + state$o.E + state$o.I)/N * 100 }, bg_color, "solid")

    p5 <- addExtraPlotQ(plots$p5, data_sample, c(dstartdate, plot_end_date),
                        function(state) { (state$y.R + state$o.R)/N * 100 }, bg_color, "solid")

    list(p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5)
}

add_all_plots <- function(plots) {
    dateRange <- c(dstartdate, plot_end_date)

    p1 <- addExtraPlotQ2(plots$p1, data_sample, c(dstartdate, plot_end_date),
                        function(state) { state$y.deadi + state$o.deadi }, "#3366FF", "solid")
    dm1 <- numeric(length(p1$data$x))
    dm1[1:length(p1$data$x)] = NA
    dm1[1:length(dmorti)] = dmorti
    p1 <- p1 + geom_line(aes(y = dm1)) + geom_point(aes(y = dm1), size=1, color="#1144CC")

    p2 <- addExtraPlotQ2(plots$p2, data_sample, c(dstartdate, plot_end_date),
                        function(state) { state$y.died + state$o.died }, "#3366FF", "solid")
    dm2 <- numeric(length(p2$data$x))
    dm2[1:length(p2$data$x)] = NA
    dm2[1:length(dmort)] = dmort
    p2 <- p2 + geom_line(aes(y = dm2)) + geom_point(aes(y = dm2),  size=1, color="#1144CC")

    p3 <- addExtraPlotQ2(plots$p3, data_sample, c(dstartdate, plot_end_date),
                        function(state) { state$y.hospi + state$o.hospi }, "#FF6633", "solid")
    dm3 <- numeric(length(p3$data$x))
    dm3[1:length(p3$data$x)] = NA
    dm3[1:length(dhospi)] = dhospi
    p3 <- p3 + geom_line(aes(y = dm3)) + geom_point(aes(y = dm3),  size=1, color="#CC4411")

    p4 <- addExtraPlotQ2(plots$p4, data_sample, c(dstartdate, plot_end_date),
                        function(state) { (state$y.E + state$y.I + state$o.E + state$o.I)/N * 100 }, "#FFFF66", "solid")

    p5 <- addExtraPlotQ2(plots$p5, data_sample, c(dstartdate, plot_end_date),
                        function(state) { (state$y.R + state$o.R)/N * 100 }, "#33FF66", "solid")

    list(p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5)
}

all_plots <- function(date_markers) {
    dateRange <- c(dstartdate, plot_end_date)

    p1 <- makePlot(data_sample, c(dstartdate, plot_end_date),
                   function(state) { state$y.deadi + state$o.deadi }, "#3366FF",
                   c("Count", "Deaths per day"), date_markers)
    dm1 <- numeric(length(p1$data$x))
    dm1[1:length(p1$data$x)] = NA
    dm1[1:length(dmorti)] = dmorti
    p1 <- p1 + geom_line(aes(y = dm1)) + geom_point(aes(y = dm1),  size=1, color="#1144CC")

    p2 <- makePlot(data_sample, c(dstartdate, plot_end_date),
                   function(state) { state$y.died + state$o.died }, "#3366FF",
                   c("Count", "Total deaths"), date_markers)

    dm2 <- numeric(length(p2$data$x))
    dm2[1:length(p2$data$x)] = NA
    dm2[1:length(dmort)] = dmort
    p2 <- p2 + geom_line(aes(y = dm2)) + geom_point(aes(y = dm2),  size=1, color="#1144CC")

    p3 <- makePlot(data_sample, c(dstartdate, plot_end_date),
                   function(state) { state$y.hospi + state$o.hospi }, "#FF6633",
                   c("Count", paste(c(HospLabel, "per day"))), date_markers)
    dm3 <- numeric(length(p3$data$x))
    dm3[1:length(p3$data$x)] = NA
    dm3[1:length(dhospi)] = dhospi
    p3 <- p3 + geom_line(aes(y = dm3)) + geom_point(aes(y = dm3),  size=1, color="#CC4411")

    p4 <- makePlot(data_sample, c(dstartdate, plot_end_date),
                   function(state) { (state$y.E + state$y.I + state$o.E + state$o.I)/N * 100 },
                   "#FFFF66", c("Belgian population (%)", "Infected people"), date_markers)

    p5 <- makePlot(data_sample, c(dstartdate, plot_end_date),
                   function(state) { (state$y.R + state$o.R)/N * 100 }, "#33FF66",
                   c("Belgian population (%)", "Recovered from disease"), date_markers)

    list(p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5)
}

all_plots_age <- function(date_markers) {
    dateRange <- c(dstartdate, plot_end_date)

    ageGroupLabels = c('< 65y','>= 65y','All (50%, 90% cri)')

    ## Plot 1 : daily deaths
    
    p1 <- makePlot(data_sample, dateRange,
                   function(state) { state$y.deadi + state$o.deadi }, "#3366FF",
                   c("Count", "Deaths per day (excl. WZC)"), date_markers)
    dm1 <- numeric(length(p1$data$x))
    dm1[1:length(p1$data$x)] = NA
    dm1[1:length(dmorti)] = dmorti
    p1 <- p1 + geom_line(aes(y = dm1)) + geom_point(aes(y = dm1),  size=1, color="#1144CC")

    ## Add y curves
    p1 <- addExtraPlot(p1, data_sample, dateRange, function(state) { state$y.deadi }, "#3366FF", "dashed")

    dm1y <- numeric(length(p1$data$x))
    dm1y[1:length(p1$data$x)] = NA
    dm1y[1:length(y.dmorti)] = y.dmorti
    
    p1 <- p1 + geom_line(aes(y = dm1y), linetype="dashed") + geom_point(aes(y = dm1y),  size=1, color="#1144CC")

    ## Add o curves
    p1 <- addExtraPlot(p1, data_sample, dateRange, function(state) { state$o.deadi }, "#3366FF", "dotted")

    dm1o <- numeric(length(p1$data$x))
    dm1o[1:length(p1$data$x)] = NA
    dm1o[1:length(o.dmorti)] = o.dmorti

    p1 <- p1 + geom_line(aes(y = dm1o), linetype="dotted") + geom_point(aes(y = dm1o),  size=1, color="#1144CC")

    p1 <- p1 + scale_colour_manual(values = c("#3366FF"), guide=FALSE) +
        scale_linetype_manual(name = 'Age group',
                              values=c('solid'='solid','dashed'='dashed','dotted'='dotted'),
                              labels = ageGroupLabels)
    
    # Plot 2 : total deaths
    
    p2 <- makePlot(data_sample, dateRange,
                   function(state) { state$y.died + state$o.died }, "#3366FF",
                   c("Count", "Total deaths (excl. WZC)"), date_markers)

    dm2 <- numeric(length(p2$data$x))
    dm2[1:length(p2$data$x)] = NA
    dm2[1:length(dmort)] = dmort
    p2 <- p2 + geom_line(aes(y = dm2)) + geom_point(aes(y = dm2), size=1, color="#1144CC")

    ## Add y curves
    p2 <- addExtraPlot(p2, data_sample, dateRange, function(state) { state$y.died }, "#3366FF", "dashed")

    dm2y <- numeric(length(p2$data$x))
    dm2y[1:length(p2$data$x)] = NA
    dm2y[1:length(y.dmort)] = y.dmort
    
    p2 <- p2 + geom_line(aes(y = dm2y), linetype="dashed") + geom_point(aes(y = dm2y),  size=1, color="#1144CC")

    ## Add o curves
    p2 <- addExtraPlot(p2, data_sample, dateRange, function(state) { state$o.died }, "#3366FF", "dotted")

    dm2o <- numeric(length(p2$data$x))
    dm2o[1:length(p2$data$x)] = NA
    dm2o[1:length(o.dmort)] = o.dmort

    p2 <- p2 + geom_line(aes(y = dm2o), linetype="dotted") + geom_point(aes(y = dm2o),  size=1, color="#1144CC")

    p2 <- p2 + scale_colour_manual(values = c("#1144CC"), guide=FALSE) +
        scale_linetype_manual(name = 'Age group',
                              values=c('solid'='solid','dashed'='dashed','dotted'='dotted'),
                              labels = ageGroupLabels)

    ## Plot 3 : hospitalisations

    p3 <- makePlot(data_sample, dateRange,
                   function(state) state$y.hospi + state$o.hospi, "#FF6633",
                   c("Count", paste(c(HospLabel, "per day"))), date_markers)
    dm3 <- numeric(length(p3$data$x))
    dm3[1:length(p3$data$x)] = NA
    dm3[1:length(dhospi)] = dhospi
    p3 <- p3 + geom_line(aes(y = dm3)) + geom_point(aes(y = dm3),  size=1, color="#CC4411")

    ## Add y/o curves
    p3 <- addExtraPlot(p3, data_sample, dateRange, function(state) { state$y.hospi }, "#FF6633", "dashed")
    p3 <- addExtraPlot(p3, data_sample, dateRange, function(state) { state$o.hospi }, "#FF6633", "dotted")

    p3 <- p3 + scale_colour_manual(values = c("#FF6633"), guide=FALSE) +
        scale_linetype_manual(name = 'Age group',
                              values=c('solid'='solid','dashed'='dashed','dotted'='dotted'),
                              labels = ageGroupLabels)

    # p3 <- p3 + scale_y_continuous(limits=c(0, 2500))

    ## Plot 4 : infected people
    
    p4 <- makePlot(data_sample, dateRange,
                   function(state) { (state$y.E + state$y.I + state$o.E + state$o.I)/N * 100 },
                   "#A67514", c("Population group (%)", "Infected people (excl. WZC)"), date_markers)

    ## Add y/o curves
    p4 <- addExtraPlot(p4, data_sample, dateRange, function(state) { (state$y.E + state$y.I)/y.N * 100 }, "#A67514", "dashed")
    p4 <- addExtraPlot(p4, data_sample, dateRange, function(state) { (state$o.E + state$o.I)/o.N * 100 }, "#A67514", "dotted")

    p4 <- p4 + scale_colour_manual(values = c("#A67514"), guide=FALSE) +
        scale_linetype_manual(name = 'Age group',
                              values=c('solid'='solid','dashed'='dashed','dotted'='dotted'),
                              labels = ageGroupLabels)
    
    ## Plot 5 : recovered people

    p5 <- makePlot(data_sample, dateRange,
                   function(state) { (state$y.R + state$o.R)/N * 100 }, "#33CC66",
                   c("Population group (%)", "Recovered from disease (excl. WZC)"), date_markers)

    ## Add y/o curves
    p5 <- addExtraPlot(p5, data_sample, dateRange, function(state) { state$y.R/y.N * 100 }, "#33CC66", "dashed")
    p5 <- addExtraPlot(p5, data_sample, dateRange, function(state) { state$o.R/o.N * 100 }, "#33CC66", "dotted")

    p5 <- p5 + scale_colour_manual(values = c("#33CC66"), guide=FALSE) +
        scale_linetype_manual(name = 'Age group',
                              values=c('solid'='solid','dashed'='dashed','dotted'='dotted'),
                              labels = ageGroupLabels)

    grid.arrange(p1, p2, p3, p4, p5, nrow=3)
}

options(scipen=999)

readSample <- function() {
    all <- readData(inputfiles)

    posterior <- all

    ## print(ess(posterior))

    ## plot(ts(subset(posterior, select=keyparamnames)))

    ## compute credibility intervals
    ##print(data.frame(ci(posterior, ci=0.01)))
    ##print(data.frame(ci(posterior, ci=0.50)))
    ##print(data.frame(ci(posterior, ci=0.95)))

    posterior

    selection <- which(posterior[,"y.IFR"] < 3)
    scount <- length(selection)
    draws <- sample(floor(scount/8):scount, quantilePlotSampleSize)
    data_sample <- posterior[selection[draws],][0:-1]

    data_sample
}

outputdir <- "output"

system(paste("mkdir ", outputdir))

quantilePlotSampleSize <- 50
data_sample <- readSample()

calcbetas.age.orig <- calcbetas.age

calcbetas.age.relax <- function(time, betay0, betayt, betao0, betaot, betayo0, betayot)
{
    if (time < relax1.start)
        return(calcbetas.age.orig(time, betay0, betayt, betao0, betaot, betayo0, betayot))

    if (time > relax3.end)
        return(c(betay0, betao0, betayo0))

    if (time > relax2.end) {
        y.beta = max(0.001,
                     relax3.measures.y.E * betayt + (1 - relax3.measures.y.E) * betay0)
        o.beta = max(0.001,
                     relax3.measures.o.E * betaot + (1 - relax3.measures.o.E) * betao0)
        yo.beta = max(0.001,
                      relax3.measures.yo.E * betayot + (1 - relax3.measures.yo.E) * betayo0)
    } else if (time > relax1.end) {
        y.beta = max(0.001,
                     relax2.measures.y.E * betayt + (1 - relax2.measures.y.E) * betay0)
        o.beta = max(0.001,
                     relax2.measures.o.E * betaot + (1 - relax2.measures.o.E) * betao0)
        yo.beta = max(0.001,
                      relax2.measures.yo.E * betayot + (1 - relax2.measures.yo.E) * betayo0)
    } else {
        y.beta = max(0.001,
                     relax1.measures.y.E * betayt + (1 - relax1.measures.y.E) * betay0)
        o.beta = max(0.001,
                     relax1.measures.o.E * betaot + (1 - relax1.measures.o.E) * betao0)
        yo.beta = max(0.001,
                      relax1.measures.yo.E * betayot + (1 - relax1.measures.yo.E) * betayo0)
    }

    c(y.beta, o.beta, yo.beta)
}

calcbetas.age <- calcbetas.age.relax

plot_end_date <- as.Date("2020/12/1")
relax1.start_date <- as.Date("2020/3/11")
relax1.end_date <- as.Date("2021/12/1")
relax2.end_date <- as.Date("2021/12/1")
relax3.end_date <- as.Date("2021/12/1")

relax1.measures.y.E = 0
relax1.measures.o.E = 0
relax1.measures.yo.E = 0

relax2.measures.y.E = 0
relax2.measures.o.E = 0
relax2.measures.yo.E = 0

relax3.measures.y.E = 0
relax3.measures.o.E = 0
relax3.measures.yo.E = 0

relax1.start <- as.numeric(relax1.start_date - dstartdate) 
relax1.end <- as.numeric(relax1.end_date - dstartdate) 
relax2.end <- as.numeric(relax2.end_date - dstartdate) 
relax3.end <- as.numeric(relax3.end_date - dstartdate)

grey95 <- "#E0E0E0"
grey75 <- "#DDDDDD"

plots <- all_plots_bg(data.frame(pos=c(as.Date("2020/03/13"), as.Date("2020/05/03")),
                                 color=c("red", "green")))

relax1.start_date <- as.Date("2021/12/1")

relax1.start <- as.numeric(relax1.start_date - dstartdate) 

plots <- add_all_plots_bg(plots)

relax1.start_date <- as.Date("2020/5/4")
relax1.end_date <- as.Date("2020/5/10")
relax2.end_date <- as.Date("2020/5/17")
relax3.end_date <- as.Date("2020/12/1")

relax1.measures.y.E = 0.95
relax1.measures.o.E = 0.95
relax1.measures.yo.E = 0.95

relax2.measures.y.E = 0.9
relax2.measures.o.E = 0.9
relax2.measures.yo.E = 0.9

relax3.measures.y.E = 0.85
relax3.measures.o.E = 0.85
relax3.measures.yo.E = 0.85

relax1.start <- as.numeric(relax1.start_date - dstartdate) 
relax1.end <- as.numeric(relax1.end_date - dstartdate) 
relax2.end <- as.numeric(relax2.end_date - dstartdate) 
relax3.end <- as.numeric(relax3.end_date - dstartdate)

plots <- add_all_plots(plots)

plots$p1 <- plots$p1 + coord_cartesian(ylim = c(0, 200), xlim=c(as.Date("2020/3/1"), as.Date("2020/7/1")))

plots$p2 <- plots$p2 + coord_cartesian(ylim = c(0, 7000), xlim=c(as.Date("2020/3/1"), as.Date("2020/7/1")))

plots$p3 <- plots$p3 + coord_cartesian(ylim = c(0, 1000), xlim=c(as.Date("2020/3/1"), as.Date("2020/7/1")))

plots$p4 <- plots$p4 + coord_cartesian(ylim = c(0, 5), xlim=c(as.Date("2020/3/1"), as.Date("2020/7/1")))

plots$p5 <- plots$p5 + coord_cartesian(ylim = c(0, 20), xlim=c(as.Date("2020/3/1"), as.Date("2020/7/1")))

svg("p1-zoom1.svg", width=8, height=6)
plots$p1
dev.off()

svg("p2-zoom1.svg", width=8, height=6)
plots$p2
dev.off()

svg("p3-zoom1.svg", width=8, height=6)
plots$p3
dev.off()

svg("p4-zoom1.svg", width=8, height=6)
plots$p4
dev.off()

svg("p5-zoom1.svg", width=8, height=6)
plots$p5
dev.off()

plots$p1 <- plots$p1 + coord_cartesian(ylim = c(0, 200), xlim=c(as.Date("2020/3/1"), as.Date("2020/12/1")))

plots$p2 <- plots$p2 + coord_cartesian(ylim = c(0, 14000), xlim=c(as.Date("2020/3/1"), as.Date("2020/12/1")))

plots$p3 <- plots$p3 + coord_cartesian(ylim = c(0, 1000), xlim=c(as.Date("2020/3/1"), as.Date("2020/12/1")))

plots$p4 <- plots$p4 + coord_cartesian(ylim = c(0, 5), xlim=c(as.Date("2020/3/1"), as.Date("2020/12/1")))

plots$p5 <- plots$p5 + coord_cartesian(ylim = c(0, 50), xlim=c(as.Date("2020/3/1"), as.Date("2020/12/1")))

svg("p1-zoom2.svg", width=8, height=6)
plots$p1
dev.off()

svg("p2-zoom2.svg", width=8, height=6)
plots$p2
dev.off()

svg("p3-zoom2.svg", width=8, height=6)
plots$p3
dev.off()

svg("p4-zoom2.svg", width=8, height=6)
plots$p4
dev.off()

svg("p5-zoom2.svg", width=8, height=6)
plots$p5
dev.off()

plots$p1 <- plots$p1 + coord_cartesian(ylim = c(0, 4000), xlim=c(as.Date("2020/3/1"), as.Date("2020/6/1")))

plots$p2 <- plots$p2 + coord_cartesian(ylim = c(0, 70000), xlim=c(as.Date("2020/3/1"), as.Date("2020/6/1")))

plots$p3 <- plots$p3 + coord_cartesian(ylim = c(0, 15000), xlim=c(as.Date("2020/3/1"), as.Date("2020/6/1")))

plots$p4 <- plots$p4 + coord_cartesian(ylim = c(0, 40), xlim=c(as.Date("2020/3/1"), as.Date("2020/6/1")))

plots$p5 <- plots$p5 + coord_cartesian(ylim = c(0, 100), xlim=c(as.Date("2020/3/1"), as.Date("2020/6/1")))

svg("p1-zoom3.svg", width=8, height=6)
plots$p1
dev.off()

svg("p2-zoom3.svg", width=8, height=6)
plots$p2
dev.off()

svg("p3-zoom3.svg", width=8, height=6)
plots$p3
dev.off()

svg("p4-zoom3.svg", width=8, height=6)
plots$p4
dev.off()

svg("p5-zoom3.svg", width=8, height=6)
plots$p5
dev.off()



#### Delaying lockdown exit

plot_end_date <- as.Date("2021/3/1")
relax1.start_date <- as.Date("2020/3/11")
relax1.end_date <- as.Date("2021/12/1")
relax2.end_date <- as.Date("2021/12/1")
relax3.end_date <- as.Date("2021/12/1")

relax1.measures.y.E = 0
relax1.measures.o.E = 0
relax1.measures.yo.E = 0

relax2.measures.y.E = 0
relax2.measures.o.E = 0
relax2.measures.yo.E = 0

relax3.measures.y.E = 0
relax3.measures.o.E = 0
relax3.measures.yo.E = 0

relax1.start <- as.numeric(relax1.start_date - dstartdate) 
relax1.end <- as.numeric(relax1.end_date - dstartdate) 
relax2.end <- as.numeric(relax2.end_date - dstartdate)  
relax3.end <- as.numeric(relax3.end_date - dstartdate)

grey95 <- "#E0E0E0"
grey75 <- "#DDDDDD"

relax1.start_date <- as.Date("2020/5/4")
relax1.end_date <- as.Date("2020/5/10")
relax2.end_date <- as.Date("2020/5/17")
relax3.end_date <- as.Date("2021/3/1")

relax1.measures.y.E = 0.95
relax1.measures.o.E = 0.95
relax1.measures.yo.E = 0.95

relax2.measures.y.E = 0.9
relax2.measures.o.E = 0.9
relax2.measures.yo.E = 0.9

relax3.measures.y.E = 0.85
relax3.measures.o.E = 0.85
relax3.measures.yo.E = 0.85

relax1.start <- as.numeric(relax1.start_date - dstartdate) 
relax1.end <- as.numeric(relax1.end_date - dstartdate) 
relax2.end <- as.numeric(relax2.end_date - dstartdate) 
relax3.end <- as.numeric(relax3.end_date - dstartdate)

plots <- all_plots_bg(data.frame(pos=c(as.Date("2020/03/13"), as.Date("2020/05/03"), as.Date("2020/05/17")),
                                 color=c("red", "green", "yellow")))

relax1.start_date <- as.Date("2020/5/17")
relax1.end_date <- as.Date("2020/5/24")
relax2.end_date <- as.Date("2020/5/31")

relax1.start <- as.numeric(relax1.start_date - dstartdate) 
relax1.end <- as.numeric(relax1.end_date - dstartdate) 
relax2.end <- as.numeric(relax2.end_date - dstartdate) 
relax3.end <- as.numeric(relax3.end_date - dstartdate)

plots <- add_all_plots(plots)
