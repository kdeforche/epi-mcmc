library(bayesplot)
library(ggplot2)
library(scales)
require(mcmcse)
require(bayestestR)
require(gridExtra)

source("control.R")
source(settings)

source(data, chdir=T)
source(fitmodel, chdir=T)

sourceR <- function(file) {
    source(paste(Rdir, file, sep=""))
}

sourceR("lib/libMCMC.R")

all_plots <- function(date_markers) {
    p1 <- makePlot(data_sample, c(dstartdate, plot_end_date),
                   function(state) state$deadi, "#3366FF",
                   c("Count", "Deaths per day"), date_markers)
    dm1 <- numeric(length(p1$data$x))
    dm1[1:length(p1$data$x)] = NA
    dm1[1:length(dmorti)] = dmorti
    p1 <- p1 + geom_line(aes(y = dm1)) + geom_point(aes(y = dm1),  size=1, color="#1144CC")

    p2 <- makePlot(data_sample, c(dstartdate, plot_end_date),
                   function(state) state$died, "#3366FF",
                   c("Count", "Total deaths"), date_markers)

    dm2 <- numeric(length(p2$data$x))
    dm2[1:length(p2$data$x)] = NA
    dm2[1:length(dmort)] = dmort
    p2 <- p2 + geom_line(aes(y = dm2)) + geom_point(aes(y = dm2),  size=1, color="#1144CC")

    p3 <- makePlot(data_sample, c(dstartdate, plot_end_date),
                   function(state) state$hospi, "#FF6633",
                   c("Count", paste(c(HospLabel, "per day"))), date_markers)
    dm3 <- numeric(length(p3$data$x))
    dm3[1:length(p3$data$x)] = NA
    dm3[1:length(dhospi)] = dhospi
    p3 <- p3 + geom_line(aes(y = dm3)) + geom_point(aes(y = dm3),  size=1, color="#CC4411")

    p4 <- makePlot(data_sample, c(dstartdate, plot_end_date),
                   function(state) { (state$E + state$I)/N * 100 },
                   "#FFFF66", c("Belgian population (%)", "Infected people"), date_markers)

    p5 <- makePlot(data_sample, c(dstartdate, plot_end_date),
                   function(state) { state$R/N * 100 }, "#33FF66",
                   c("Belgian population (%)", "Recovered from disease"), date_markers)

    grid.arrange(p1, p2, p3, p4, p5, nrow=3)
}

all_plots_age <- function(date_markers) {
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
                   function(state) state$y.hospi + state$o.hospi, "#FF6633",
                   c("Count", paste(c(HospLabel, "per day"))), date_markers)
    dm3 <- numeric(length(p3$data$x))
    dm3[1:length(p3$data$x)] = NA
    dm3[1:length(dhospi)] = dhospi
    p3 <- p3 + geom_line(aes(y = dm3)) + geom_point(aes(y = dm3),  size=1, color="#CC4411")

    p4 <- makePlot(data_sample, c(dstartdate, plot_end_date),
                   function(state) { (state$y.E + state$y.I + state$o.E + state$y.I)/N * 100 },
                   "#FFFF66", c("Belgian population (%)", "Infected people"), date_markers)

    p5 <- makePlot(data_sample, c(dstartdate, plot_end_date),
                   function(state) { (state$y.R + state$o.R)/N * 100 }, "#33FF66",
                   c("Belgian population (%)", "Recovered from disease"), date_markers)

    grid.arrange(p1, p2, p3, p4, p5, nrow=3)
}

all_plots_jo <- function(date_markers) {
    p1 <- makePlot(data_sample, c(dstartdate, plot_end_date),
                   function(state) { state$j.deadi + state$o.deadi }, "#3366FF",
                   c("Count", "Deaths per day"), date_markers)

    p2 <- makePlot(data_sample, c(dstartdate, plot_end_date),
                   function(state) { state$j.died + state$o.died }, "#3366FF",
                   c("Count", "Total deaths"), date_markers)

    p3 <- makePlot(data_sample, c(dstartdate, plot_end_date),
                   function(state) state$j.hospi + state$o.hospi, "#FF6633",
                   c("Count", paste(c(HospLabel, "per day"))), date_markers)

    p4 <- makePlot(data_sample, c(dstartdate, plot_end_date),
                   function(state) { (state$j.E + state$j.I + state$o.E + state$o.I)/N * 100 },
                   "#FFFF66", c("Belgian population (%)", "Infected people"), date_markers)

    p5 <- makePlot(data_sample, c(dstartdate, plot_end_date),
                   function(state) { (state$j.R + state$o.R)/N * 100 }, "#33FF66",
                   c("Belgian population (%)", "Recovered from disease"), date_markers)

    grid.arrange(p1, p2, p3, p4, p5, nrow=3)
}

options(scipen=999)

system(paste("mkdir ", outputdir))

# read new data
all <- readData(inputfiles)

posterior <- all

# calculate effective sample sizes.
pess <- ess(posterior)
pess

plot(ts(subset(posterior, select=c("y.R0", "y.Rt", "IFR", "Tinf", "betay0", "betayt"))))

# compute credibility intervals
ci(posterior, method = "HDI", ci=0.01)
ci(posterior, method = "HDI", ci=0.50)
ci(posterior, method = "HDI", ci=0.95)

###### Create sample

## sample a subset of data to show in density plots
selection <- which(posterior[,"IFR"] < 3)
scount <- length(selection)
draws <- sample(floor(scount/8):scount, quantilePlotSampleSize)
data_sample <- posterior[selection[draws],][0:-1]

###### Current state

all_plots <- all_plots_age

##sourceR("models/model.R")

pdf(paste(outputdir, "/current-state.pdf", sep=""), width=12, height=16)

plot_end_date <- as.Date("2020/8/1")
all_plots(data.frame(pos=c(Sys.Date()), color=c("red")))

dev.off()

png(paste(outputdir, "/forecast-time.png", sep=""), width=1200, height=1600)

plot_end_date <- as.Date("2020/5/1")
all_plots(data.frame(pos=c(Sys.Date()), color=c("red")))

dev.off()

######## Simple exiting scenario's

sourceR("models/model-exiting-lockdown.R")

plot_end_date <- as.Date("2020/11/1")

labels <- c("0", "01", "02", "03", "04", "05", "06", "07", "08", "09")

i = 1
for (r in seq(0,0.9,0.1)) {
    relax_date <- as.Date("2020/6/30")
    lift_date <- as.Date("2020/8/31")

    relax_E = r # 0 : no lockdown, 1 : lockdown 'lite' as now
    relax_offset <- as.numeric(relax_date - dstartdate) 
    lift_offset <- as.numeric(lift_date - dstartdate) 

    if (relax_E == 0)
        lift_date = relax_date
    
    pdf(paste(outputdir, "/simple-relax-exit", labels[i], "-", lift_offset, ".pdf", sep=""),
        width=12, height=16)

    all_plots(data.frame(pos=c(relax_date, lift_date), color=c("#888800", "#008800")))

    i = i + 1
    
    dev.off()
}

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

######## compartiment model
quit() # doesn't work yet

sourceR("models/model-compartiment.R")

## params <- data_sample[1,]
## params <- transformParams(unlist(params, use.names=FALSE))
## stateJO <- calculateModel(Initial, params, 300)
## plot(stateJO$o.died, type='l', col='red')
## lines(stateJO$j.died, type='l', col='blue')

plot_end_date <- as.Date("2020/11/1")

## date at which we start with this
start_date <- as.Date("2020/4/30")
start_offset <- as.numeric(start_date - dstartdate)

## date at which we go out of lockdown
end_date <- as.Date("2020/8/31")
end_offset <- as.numeric(end_date - dstartdate)

o.beta0.factor <- 0.5
j.beta0.factor <- (N - o.beta0.factor * o.N) / j.N

o.betat.factor <- 0.7
j.betat.factor <- (N - o.betat.factor * o.N) / j.N

o.E <- 1.0
j.E <- 0.7
jo.E <- 1.0

pdf(paste(outputdir, "/age-compart-05-07-1-02-1.pdf", sep=""), width=12, height=16)
all_plots_jo(data.frame(pos=c(start_date, end_date), color=c("#888800", "#008800")))
dev.off()

## plots per compartement, not used yet.
##  deaths per day
## makePlot(data_sample, c(dstartdate, plot_end_date), function(state) state$j.deadi, "#3366FF", c("Count", "Deaths per day"))

## makePlot(data_sample, c(dstartdate, plot_end_date), function(state) state$o.deadi, "#3366FF", c("Count", "Deaths per day"))

## makePlot(data_sample, c(dstartdate, plot_end_date), function(state) (state$j.deadi + state$o.deadi), "#3366FF", c("Count", "Deaths per day"))

## ##  cumulative deaths
## makePlot(data_sample, c(dstartdate, plot_end_date), function(state) state$j.died, "#3366FF", c("Count", "Total deaths"))

## makePlot(data_sample, c(dstartdate, plot_end_date), function(state) state$o.died, "#3366FF", c("Count", "Total deaths"))

## makePlot(data_sample, c(dstartdate, plot_end_date), function(state) { state$j.died + state$o.died }, "#3366FF", c("Count", "Total deaths"))

## ##  hospitalisations
## makePlot(data_sample, c(dstartdate, plot_end_date), function(state) state$j.hospi, "#FF6633",  c("Count", "New hospitalisations per day"))

## makePlot(data_sample, c(dstartdate, plot_end_date), function(state) state$o.hospi, "#FF6633",  c("Count", "New hospitalisations per day"))

## makePlot(data_sample, c(dstartdate, plot_end_date), function(state) (state$j.hospi + state$o.hospi), "#FF6633",  c("Count", "New hospitalisations per day"))

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


