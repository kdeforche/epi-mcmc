library(bayesplot)
library(ggplot2)
library(scales)
require(mcmcse)
require(bayestestR)
require(gridExtra)

source("settings.R")

source(data, chdir=T)

trimData <- function(count) {
    y.dmorti <<- y.dmorti[1:(length(y.dmorti) - count)]
    o.dmorti <<- o.dmorti[1:(length(o.dmorti) - count)]
    dmorti <<- dmorti[1:(length(dmorti) - count)]

    y.dhospi <<- y.dhospi[1:(length(y.dhospi) - count)]
    o.dhospi <<- o.dhospi[1:(length(o.dhospi) - count)]
    dhospi <<- dhospi[1:(length(dhospi) - count)]

    y.dmort <<- y.dmort[1:(length(y.dmort) - count)]
    o.dmort <<- o.dmort[1:(length(o.dmort) - count)]
    dmort <<- dmort[1:(length(dmort) - count)]

    dhosp <<- dhosp[1:(length(dhosp) - count)]
}

source(fitmodel, chdir=T)

sourceR <- function(file) {
    source(paste(Rdir, file, sep=""))
}

sourceR("lib/libMCMC.R")

all_plots <- function(date_markers) {
    p1 <- makePlot(data_sample, c(dstartdate, plot_end_date),
                   function(state) state$deadi, "#3366FF",
                   c("Count", "Deaths per day"), date_markers, NULL)
    dm1 <- numeric(length(p1$data$x))
    dm1[1:length(p1$data$x)] = NA
    dm1[1:length(dmort)] = dmorti
    p1 <- p1 + geom_line(aes(y = dm1)) + geom_point(aes(y = dm1), size=0.5, color="#1144CC")
    
    p2 <- makePlot(data_sample, c(dstartdate, plot_end_date),
                   function(state) state$died, "#3366FF",
                   c("Count", "Total deaths"), date_markers, NULL)

    dm2 <- numeric(length(p2$data$x))
    dm2[1:length(p2$data$x)] = NA
    dm2[1:length(dmort)] = dmort
    p2 <- p2 + geom_line(aes(y = dm2)) + geom_point(aes(y = dm2),  size=0.5, color="#1144CC")

    p3 <- makePlot(data_sample, c(dstartdate, plot_end_date),
                   function(state) state$hospi, "#FF6633",
                   c("Count", paste(c(HospLabel, "per day"))), date_markers, NULL)
    dm3 <- numeric(length(p3$data$x))
    dm3[1:length(p3$data$x)] = NA
    dm3[1:length(dhospi)] = dhospi
    p3 <- p3 + geom_line(aes(y = dm3)) + geom_point(aes(y = dm3),  size=0.5, color="#CC4411")

    p4 <- makePlot(data_sample, c(dstartdate, plot_end_date),
                   function(state) { (state$E + state$In + state$Is)/N * 100 },
                   "#FFFF66", c(paste(country_adjective, "population (%)"), "Infected people"), date_markers, NULL)

    ## p4 <- makePlot(data_sample, c(dstartdate, plot_end_date),
    ##                function(state) { (state$E + state$I)/N * 100 },
    ##                "#FFFF66", c(paste(country_adjective, "population (%)"), "Infected people"), date_markers, NULL)

    p5 <- makePlot(data_sample, c(dstartdate, plot_end_date),
                   function(state) { state$R/N * 100 }, "#33FF66",
                   c(paste(country_adjective, "population (%)"), "Recovered from disease"), date_markers, NULL)

    p6 <- makePlot(data_sample, c(dstartdate, plot_end_date),
                   function(state) { state$Re }, "#33FFFF",
                   c("Re", "Effective reproduction number (Re)"), date_markers, NULL)

    p6 <- p6 + coord_cartesian(ylim = c(0, NA)) +
        geom_hline(yintercept=1, linetype="solid", color="gray", size=0.5)

    ## p7 <- makePlot(data_sample, c(dstartdate, plot_end_date),
    ##                function(state) { state$BS }, "#33FF66",
    ##                c(paste(country_adjective, "BS"), ""), date_markers, NULL)

    ## p8 <- makePlot(data_sample, c(dstartdate, plot_end_date),
    ##                function(state) { state$Ib * state$BS }, "#33FF66",
    ##                c(paste(country_adjective, "Ib * BS"), ""), date_markers, NULL)

    ## p9 <- makePlot(data_sample, c(dstartdate, plot_end_date),
    ##                function(state) { (state$Ib / (N / state$BS)) }, "#33FF66",
    ##                c("homogeneity", ""), date_markers, NULL)
   
    ## grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow=3)
    grid.arrange(p1, p2, p3, p4, p5, p6, nrow=3)
}

all_plots_age <- function(date_markers) {
    dateRange <- c(plot_start_date, plot_end_date)

    ageGroupLabels = c('< 65y','>= 65y','All (50%, 90% cri)')

    ## Plot 1 : daily deaths

    p1 <- makePlot(data_sample, dateRange,
                   function(state) { state$y.deadi + state$o.deadi }, "#3366FF",
                   c("Count", "Deaths per day"), date_markers, 'solid')

    start <- as.numeric(dstartdate - dateRange[1])
   
    dm1 <- numeric(length(p1$data$x))
    dm1[1:length(p1$data$x)] = NA
    dm1[(start + 1):(start + length(dmorti))] = dmorti
    p1 <- p1 + geom_line(aes(y = dm1)) + geom_point(aes(y = dm1),  size=0.5, color="#1144CC")

    if (zoom) {
        p1 <- p1 + coord_cartesian(ylim = c(0, 50))
    }

    ## Add y curves
    p1 <- addExtraPlot(p1, data_sample, dateRange, function(state) { state$y.deadi }, "#3366FF", "dashed")

    dm1y <- numeric(length(p1$data$x))
    dm1y[1:length(p1$data$x)] = NA
    dm1y[(start + 1):(start + length(y.dmorti))] = y.dmorti
    
    p1 <- p1 + geom_line(aes(y = dm1y), linetype="dashed") + geom_point(aes(y = dm1y),  size=0.5, color="#1144CC")

    ## Add o curves
    p1 <- addExtraPlot(p1, data_sample, dateRange, function(state) { state$o.deadi }, "#3366FF", "dotted")

    dm1o <- numeric(length(p1$data$x))
    dm1o[1:length(p1$data$x)] = NA
    dm1o[(start + 1):(start + length(o.dmorti))] = o.dmorti

    p1 <- p1 + geom_line(aes(y = dm1o), linetype="dotted") + geom_point(aes(y = dm1o),  size=0.5, color="#1144CC")

    p1 <- p1 + scale_colour_manual(values = c("#3366FF"), guide=FALSE) +
        scale_linetype_manual(name = 'Age group',
                              values=c('solid'='solid','dashed'='dashed','dotted'='dotted'),
                              labels = ageGroupLabels)
    
    # Plot 2 : total deaths
    
    p2 <- makePlot(data_sample, dateRange,
                   function(state) { state$y.died + state$o.died }, "#3366FF",
                   c("Count", "Total deaths"), date_markers, 'solid')

    dm2 <- numeric(length(p2$data$x))
    dm2[1:length(p2$data$x)] = NA
    dm2[(start + 1):(start + length(dmort))] = dmort
    p2 <- p2 + geom_line(aes(y = dm2)) + geom_point(aes(y = dm2), size=0.5, color="#1144CC")

    ## Add y curves
    p2 <- addExtraPlot(p2, data_sample, dateRange, function(state) { state$y.died }, "#3366FF", "dashed")

    dm2y <- numeric(length(p2$data$x))
    dm2y[1:length(p2$data$x)] = NA
    dm2y[(start + 1):(start + length(y.dmort))] = y.dmort
    
    p2 <- p2 + geom_line(aes(y = dm2y), linetype="dashed") + geom_point(aes(y = dm2y),  size=0.5, color="#1144CC")

    ## Add o curves
    p2 <- addExtraPlot(p2, data_sample, dateRange, function(state) { state$o.died }, "#3366FF", "dotted")

    dm2o <- numeric(length(p2$data$x))
    dm2o[1:length(p2$data$x)] = NA
    dm2o[(start + 1):(start + length(o.dmort))] = o.dmort

    p2 <- p2 + geom_line(aes(y = dm2o), linetype="dotted") + geom_point(aes(y = dm2o),  size=0.5, color="#1144CC")

    p2 <- p2 + scale_colour_manual(values = c("#1144CC"), guide=FALSE) +
        scale_linetype_manual(name = 'Age group',
                              values=c('solid'='solid','dashed'='dashed','dotted'='dotted'),
                              labels = ageGroupLabels)

    ## Plot 3 : hospitalisations

    p3 <- makePlot(data_sample, dateRange,
                   function(state) state$y.hospi + state$o.hospi, "#FF6633",
                   c("Count", paste(c(HospLabel, "per day"))), date_markers, 'solid')
    dm3 <- numeric(length(p3$data$x))
    dm3[1:length(p3$data$x)] = NA
    dm3[(start + 1):(start + length(dhospi))] = dhospi
    p3 <- p3 + geom_line(aes(y = dm3)) + geom_point(aes(y = dm3),  size=0.5, color="#CC4411")

    if (zoom) {
        p3 <- p3 + coord_cartesian(ylim = c(0, 800))
    }

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
                   "#A67514", c("Population group (%)", "Infected people"), date_markers, 'solid')

    if (zoom) {
        p4 <- p4 + coord_cartesian(ylim = c(0, 0.5))
    }

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
                   c("Population group (%)", "Recovered from disease"), date_markers, 'solid')

    ## Add y/o curves
    p5 <- addExtraPlot(p5, data_sample, dateRange, function(state) { state$y.R/y.N * 100 }, "#33CC66", "dashed")
    p5 <- addExtraPlot(p5, data_sample, dateRange, function(state) { state$o.R/o.N * 100 }, "#33CC66", "dotted")

    p5 <- p5 + scale_colour_manual(values = c("#33CC66"), guide=FALSE) +
        scale_linetype_manual(name = 'Age group',
                              values=c('solid'='solid','dashed'='dashed','dotted'='dotted'),
                              labels = ageGroupLabels)

    p6 <- makePlot(data_sample, c(dstartdate, plot_end_date),
                   function(state) { state$Re }, "#33FFFF",
                   c("Re", "Effective reproduction number (Re)"), date_markers, NULL)

    if (zoom) {
        p6 <- p6 + coord_cartesian(ylim = c(0, 2))
    }

    p6 <- p6 +  geom_hline(yintercept=1, linetype="solid", color="gray", size=0.5)

    
    ##grid.arrange(p1, p2, p3, p4, p5, p6, nrow=3)
    grid.arrange(p1, p3, p4, p6, nrow=2)
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

    posterior <- all[0:-7]

    plot(ts(subset(posterior, select=keyparamnames)))

    ## compute credibility intervals
    print(data.frame(ci(posterior, ci=0.01)))
    print(data.frame(ci(posterior, ci=0.50)))
    print(data.frame(ci(posterior, ci=0.95)))

    posterior

    selection <- 1:dim(posterior)[1]
    scount <- length(selection)
    draws <- sample(floor(scount/8):scount, quantilePlotSampleSize)
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

zoom <- T
pdf("current-state-zoom.pdf", width=12, height=12)

zoom <- F
pdf("current-state.pdf", width=12, height=12)

plot_start_date <- as.Date("2020/2/15")
plot_end_date <- as.Date("2020/9/1")

dates <- data.frame(pos=c(dstartdate + lockdown_offset,
                          dstartdate + lockdown_offset + lockdown_transition_period,
                          dstartdate + d3,
                          dstartdate + d4),
                    color=c("orange", "red", "green", "green"))

all_plots(dates)

## current Re value
est.Re <- data.frame(quantileData(data_sample, function(state) { state$Re }, lockdown_offset, lockdown_offset + 200, c(0.05, 0.5, 0.95)))
colnames(est.Re) <- c("q5", "q50", "q95")

now.est.Re.median <- est.Re$q50[Sys.Date() - dstartdate]
now.est.Re.cri95lo <- est.Re$q5[Sys.Date() - dstartdate]
now.est.Re.cri95hi <- est.Re$q95[Sys.Date() - dstartdate]

dev.off()
