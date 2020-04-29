library(bayesplot)
library(ggplot2)
library(scales)
require(mcmcse)
require(bayestestR)
require(gridExtra)

source("control.R")
source(settings)

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
    dm1[1:length(dmorti)] = dmorti
    p1 <- p1 + geom_line(aes(y = dm1)) + geom_point(aes(y = dm1),  size=1, color="#1144CC")

    p2 <- makePlot(data_sample, c(dstartdate, plot_end_date),
                   function(state) state$died, "#3366FF",
                   c("Count", "Total deaths"), date_markers, NULL)

    dm2 <- numeric(length(p2$data$x))
    dm2[1:length(p2$data$x)] = NA
    dm2[1:length(dmort)] = dmort
    p2 <- p2 + geom_line(aes(y = dm2)) + geom_point(aes(y = dm2),  size=1, color="#1144CC")

    p3 <- makePlot(data_sample, c(dstartdate, plot_end_date),
                   function(state) state$hospi, "#FF6633",
                   c("Count", paste(c(HospLabel, "per day"))), date_markers, NULL)
    dm3 <- numeric(length(p3$data$x))
    dm3[1:length(p3$data$x)] = NA
    dm3[1:length(dhospi)] = dhospi
    p3 <- p3 + geom_line(aes(y = dm3)) + geom_point(aes(y = dm3),  size=1, color="#CC4411")

    p4 <- makePlot(data_sample, c(dstartdate, plot_end_date),
                   function(state) { (state$E + state$I)/N * 100 },
                   "#FFFF66", c(paste(country_adjective, "population (%)"), "Infected people"), date_markers, NULL)

    p5 <- makePlot(data_sample, c(dstartdate, plot_end_date),
                   function(state) { state$R/N * 100 }, "#33FF66",
                   c(paste(country_adjective, "population (%)"), "Recovered from disease"), date_markers, NULL)

    grid.arrange(p1, p2, p3, p4, p5, nrow=3)
}

all_plots_age <- function(date_markers) {
    dateRange <- c(dstartdate, plot_end_date)

    ageGroupLabels = c('< 65y','>= 65y','All (50%, 90% cri)')

    ## Plot 1 : daily deaths
    
    p1 <- makePlot(data_sample, dateRange,
                   function(state) { state$y.deadi + state$o.deadi }, "#3366FF",
                   c("Count", "Deaths per day (excl. WZC)"), date_markers, 'solid')
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
                   c("Count", "Total deaths (excl. WZC)"), date_markers, 'solid')

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
                   c("Count", paste(c(HospLabel, "per day"))), date_markers, 'solid')
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
                   "#A67514", c("Population group (%)", "Infected people (excl. WZC)"), date_markers, 'solid')

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
                   c("Population group (%)", "Recovered from disease (excl. WZC)"), date_markers, 'solid')

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

    print(ess(posterior))

    plot(ts(subset(posterior, select=keyparamnames)))

    ## compute credibility intervals
    print(data.frame(ci(posterior, ci=0.01)))
    print(data.frame(ci(posterior, ci=0.50)))
    print(data.frame(ci(posterior, ci=0.95)))

    posterior

    selection <- 1:dim(posterior)[1]
    scount <- length(selection)
    draws <- sample(floor(scount/8):scount, quantilePlotSampleSize)
    data_sample <- posterior[selection[draws],][0:-1]

    data_sample
}

outputdir <- "output"

system(paste("mkdir ", outputdir))

# trimData(4)

quantilePlotSampleSize <- 1500
data_sample <- readSample()

## Configure this depending on the model
all_plots <- all_plots_age

pdf(paste(outputdir, "/current-state.pdf", sep=""), width=12, height=16)

plot_end_date <- as.Date("2020/8/1")
all_plots(data.frame(pos=c(as.Date("2020/4/28")), color=c("red")))

dev.off()

png(paste(outputdir, "/forecast-time.png", sep=""), width=1200, height=1600)

all_plots(data.frame(pos=c(Sys.Date()), color=c("red")))

dev.off()

######## No measures were taken

sourceR("models/model-exiting-age-groups.R")

calcbetas.age <- calcbetas.age.relax

plot_end_date <- as.Date("2020/6/1")
relax.start_date <- as.Date("2020/3/11")
relax.end_date <- as.Date("2020/12/1")

relax.measures.y.E = 0
relax.measures.o.E = 0
relax.measures.yo.E = 0
relax.start <- as.numeric(relax.start_date - dstartdate) 
relax.end <- as.numeric(relax.end_date - dstartdate) 

pdf(paste(outputdir, "/no-lockdown.pdf", sep=""), width=12, height=16)
all_plots(data.frame(pos=c(Sys.Date()), color=c("red")))
dev.off()

######## Simple exiting scenario's

sourceR("models/model-exiting-age-groups.R")

plot_end_date <- as.Date("2020/11/1")

## quantilePlotSampleSize <- 50

calcbetas.age <- calcbetas.age.relax

labels <- c("0", "01", "02", "03", "04", "05", "06", "07", "08", "09")

r <- 0.85
r <- 0.8
r <- 0

i = 1
for (r in seq(0,0.9,0.1)) {
    # relax.start_date <- as.Date("2020/4/16")
    relax.start_date <- as.Date("2020/5/4")
    relax.end_date <- as.Date("2020/12/1")

    relax.measures.y.E = r
    relax.start <- as.numeric(relax.start_date - dstartdate) 
    relax.end <- as.numeric(relax.end_date - dstartdate) 

    ##pdf(paste(outputdir, "/full-exit-zoom.pdf", sep=""), width=12, height=16)
    pdf(paste(outputdir, "/young-relax0-exit.pdf", sep=""), width=12, height=16)
    ##pdf(paste(outputdir, "/young-relax-exit", labels[i], "-",
    ##          relax.start, "-", relax.end, ".pdf", sep=""), width=12, height=16)

    all_plots(data.frame(pos=c(relax.start_date, relax.end_date), color=c("#888800", "#008800")))

    dev.off()

    i = i + 1
}

######## Experiment scenario's

## One weekend relaxed

plot_end_date <- as.Date("2020/6/1")

## quantilePlotSampleSize <- 50

calcbetas.age.experiment <- function(time, betay0, betayt, betao0, betaot, betayo0, betayot)
{
    if (time < relax.start)
        return(calcbetas.age.orig(time, betay0, betayt, betao0, betaot, betayo0, betayot))

    if (time > relax.end)
        return(calcbetas.age.orig(time, betay0, betayt, betao0, betaot, betayo0, betayot))
    
    y.beta = max(0.001, relax.measures.y.E * betayt + (1 - relax.measures.y.E) * betay0)
    o.beta = max(0.001, relax.measures.o.E * betaot + (1 - relax.measures.o.E) * betao0)
    yo.beta = max(0.001, relax.measures.yo.E * betayot + (1 - relax.measures.yo.E)
                  * betayo0)

    c(y.beta, o.beta, yo.beta)
}


calcbetas.age <- calcbetas.age.experiment

relax.start_date <- as.Date("2020/5/3")
relax.end_date <- as.Date("2020/5/5")

relax.measures.y.E = 0
relax.measures.o.E = 0
relax.measures.yo.E = 0

relax.start <- as.numeric(relax.start_date - dstartdate) 
relax.end <- as.numeric(relax.end_date - dstartdate) 

pdf(paste(outputdir, "/experiment", sep=""), width=12, height=16)

all_plots(data.frame(pos=c(relax.start_date, relax.end_date), color=c("#888800", "#008800")))
    
dev.off()



######## Lockdowns Wuhan style

## (1) 1.1

plot_end_date <- as.Date("2020/6/10")
relax_date <- as.Date("2020/5/1")
lift_date <- as.Date("2020/5/30")

relax_E = 1.1  # 0 : no lockdown, 1 : lockdown 'lite' as now
relax_offset <- as.numeric(relax_date - dstartdate) 
lift_offset <- as.numeric(lift_date - dstartdate) 

if (relax_E == 0)
    lift_date = relax_date

pdf(paste(outputdir, "/wuhan11.pdf", sep=""), width=12, height=16) 
all_plots(data.frame(pos=c(relax_date, lift_date), color=c("#888800", "#008800")))
dev.off()

## (2) 1.05

plot_end_date <- as.Date("2020/6/10")
relax_date <- as.Date("2020/5/1")
lift_date <- as.Date("2020/5/30")

relax_E = 1.05  # 0 : no lockdown, 1 : lockdown 'lite' as now
relax_offset <- as.numeric(relax_date - dstartdate) 
lift_offset <- as.numeric(lift_date - dstartdate) 

if (relax_E == 0)
    lift_date = relax_date

pdf(paste(outputdir, "/wuhan105.pdf", sep=""), width=12, height=16) 
all_plots(data.frame(pos=c(relax_date, lift_date), color=c("#888800", "#008800")))
dev.off()

##
## Spain exiting scenario
##
sourceR("models/model.R")
sourceR("models/model-exiting-lockdown.R")

plot_end_date <- as.Date("2020/11/1")

relax_date <- as.Date("2020/4/14")
lift_date <- as.Date("2020/8/31")

relax_E = 0.6 # 0 : no lockdown, 1 : lockdown 'lite' as now
relax_offset <- as.numeric(relax_date - dstartdate) 
lift_offset <- as.numeric(lift_date - dstartdate) 

if (relax_E == 0)
    lift_date = relax_date

pdf(paste(outputdir, "/spain-exit-0.6.pdf", sep=""), width=12, height=16)

all_plots(data.frame(pos=c(relax_date, lift_date), color=c("#888800", "#008800")))

dev.off()

png(paste(outputdir, "/spain-exit-0.6.png", sep=""), width=1200, height=1600)

all_plots(data.frame(pos=c(relax_date, lift_date), color=c("#888800", "#008800")))

dev.off()


###
