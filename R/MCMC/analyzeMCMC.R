library(bayesplot)
library(ggplot2)
library(scales)
require(mcmcse)
require(bayestestR)
require(gridExtra)

## Read settings file
source("settings.R")

## Load data, and fitmodel (defined in control.R)
source(data, chdir=T)
source(fitmodel, chdir=T)

sourceR <- function(file) {
    source(paste(Rdir, file, sep=""))
}

sourceR("lib/libMCMC.R")
sourceR("lib/libfit.R")

source("settings.R")

options(scipen=999)

DIC <- function(posterior1) {
  draw = posterior1

  x = draw[0:-7]
  lik = draw$loglikelihood
  lik.fun = function(x) { calclogl(transformParams(x)) }

  dmorti <<- dmorti[1:(length(dmorti) - evaluation_data_count)]
  dmort <<- dmort[1:(length(dmorti) - evaluation_data_count)]
  dhospi <<- dhospi[1:(length(dhospi) - evaluation_data_count)]
  dhosp <<- dhosp[1:(length(dhosp) - evaluation_data_count)]
  
  D.bar <- -2*mean(lik)
  if(is.vector(x)) theta.bar = mean(x) else theta.bar <- apply(x,2,mean)
  D.hat <- -2*lik.fun(theta.bar)
  pD <- D.bar - D.hat
  pV <- var(-2*lik)/2
  list(DIC=pD+D.bar,DIC2=pV+D.bar,IC=2*pD+D.bar,pD=pD,pV=pV,Dbar=D.bar,Dhat=D.hat)
}

evaluation_data_count <- max(0, (dstartdate + length(dmorti)) - (d2 + fitPeriod))
print(evaluation_data_count)

all_plots <- function(date_markers, range) {
    p1 <- makePlot(data_sample, range,
                   function(state, params) state$deadi, "#3366FF",
                   c("Death incidence", "Incidence of deaths"), date_markers, NULL)

    dm1c = length(dmorti) - evaluation_data_count

    start <- as.numeric(dstartdate - range[1])
    
    dm1 <- numeric(length(p1$data$x))
    dm1[1:length(p1$data$x)] = NA
    dm1[(start + 1):(start + dm1c)] = dmorti[1:dm1c]
    p1 <- p1 + geom_line(aes(y = dm1)) + geom_point(aes(y = dm1),  size=1, color="#1144CC")

    if (evaluation_data_count > 0) {
        dm1e <- numeric(length(p1$data$x))
        dm1e[1:length(p1$data$x)] = NA
        dm1e[(start + dm1c + 1):(start + length(dmorti))] = dmorti[(dm1c + 1):length(dmorti)]
        p1 <- p1 + geom_line(aes(y = dm1e), linetype='dashed', color="#555555") +
            geom_point(aes(y = dm1e), size=1, color="#555555")
    }

    p2 <- makePlot(data_sample, range,
                   function(state, params) state$died, "#3366FF",
                   c("Count", "Total deaths"), date_markers, NULL)

    dm2c = length(dmort) - evaluation_data_count

    dm2 <- numeric(length(p2$data$x))
    dm2[1:length(p2$data$x)] = NA
    dm2[(start + 1):(start + dm2c)] = dmort[1:dm2c]
    p2 <- p2 + geom_line(aes(y = dm2)) + geom_point(aes(y = dm2),  size=1, color="#1144CC")

    if (evaluation_data_count) {
        dm2e <- numeric(length(p2$data$x))
        dm2e[1:length(p2$data$x)] = NA
        dm2e[(start + dm2c + 1):(start + length(dmort))] = dmort[(dm2c + 1):length(dmort)]
        p2 <- p2 + geom_line(aes(y = dm2e), linetype='dashed', color="#555555") +
            geom_point(aes(y = dm2e), size=1, color="#555555")
    }

    p3 <- makePlot(data_sample, range,
                   function(state, params) state$hospi, "#FF6633",
                   c("Count", paste(c(HospLabel, "per day"))), date_markers, NULL)

    dm3c = length(dhospi) - evaluation_data_count

    dm3 <- numeric(length(p3$data$x))
    dm3[1:length(p3$data$x)] = NA
    dm3[(start + 1):(start + dm3c)] = dhospi[1:dm3c]
    p3 <- p3 + geom_line(aes(y = dm3)) + geom_point(aes(y = dm3),  size=1, color="#CC4411")

    if (evaluation_data_count > 0) {
        dm3e <- numeric(length(p3$data$x))
        dm3e[1:length(p3$data$x)] = NA
        dm3e[(start + dm3c + 1):(start + length(dhospi))] = dhospi[(dm3c + 1):length(dhospi)]
        p3 <- p3 + geom_line(aes(y = dm3e), linetype='dashed', color="#555555") +
            geom_point(aes(y = dm3e), size=1, color="#555555")
    }

    p4 <- makePlot(data_sample, range,
                   function(state, params) { if ("In" %in% names(state)) {
                                         return ((state$E + state$In + state$Is)/N * 100)
                                     } else {
                                         return ((state$E + state$I)/N * 100)
                                     }
                   }, "#FFFF66", c(paste("% of population of", country_adjective), "Infected individuals (%)"), date_markers, NULL)

    p5 <- makePlot(data_sample, range,
                   function(state, params) { state$R/N * 100 }, "#33FF66",
                   c(paste("% of population of", country_adjective), "Removed (%)"), date_markers, NULL)

    p6 <- makePlot(data_sample, range,
                   function(state, params) { state$Re }, "#33FFFF",
                   c("Re", "Effective reproduction number (Re)"), date_markers, NULL)

    p6 <- p6 + coord_cartesian(ylim = c(0, NA)) +
        geom_hline(yintercept=1, linetype="solid", color="#11AA11", size=0.5)

    ## Plot all mobility data
    
    
    grid.arrange(p1, p2, p3, p4, p5, p6, nrow=3)
}

read <- function() {
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

    posterior <- all

    ##print(ess(posterior))

    pdf(paste(country2,"sample.pdf",sep='_'), width=12, height=16)
    plot(ts(subset(posterior, select=keyparamnames)))
    dev.off()

    ## compute credibility intervals
    print(data.frame(ci(posterior, ci=0.01)))
    print(data.frame(ci(posterior, ci=0.50)))
    print(data.frame(ci(posterior, ci=0.95)))

    posterior
}

posterior1 <- read()
posterior <- posterior1[0:-7]

takeSample <- function(posterior) {
    ##print(ess(posterior))

    ## compute credibility intervals
    selection <- 1:dim(posterior)[1]
    scount <- length(selection)
    draws <- sample(floor(scount/8):scount, quantilePlotSampleSize)
    data_sample <- posterior[selection[draws],]

    data_sample
}

data_sample <- takeSample(posterior)

pdf(paste(country2,"graphs.pdf",sep='_'), width=12, height=16)
plot_start_date <- min(as.Date("2020/2/15"), dstartdate)
plot_end_date <- as.Date("2020/7/1")
##plot_end_date <- as.Date("2020/11/1")
##plot_end_date <- as.Date("2020/7/15")
range <- c(plot_start_date, plot_end_date)
## dates <- data.frame(pos=c(dstartdate + lockdown_offset,
##                           dstartdate + lockdown_offset + lockdown_transition_period,
##                           dstartdate + relax_offset,
##                           dstartdate + lockdown2_offset),
##                     color=c("orange", "red", "green", "red"))
dates <- data.frame(pos=c(dstartdate + lockdown_offset,
                          dstartdate + lockdown_offset + lockdown_transition_period),
                    color=c("orange", "red"))
all_plots(dates, range)
dev.off()

result <- data.frame(country = c(country2))
result$N <- N
result$lockdown.d1 <- dstartdate + lockdown_offset
result$lockdown.d2 <- dstartdate + lockdown_offset + lockdown_transition_period
result$d1.deaths <- dmort[lockdown_offset]
result$d2.deaths <- dmort[lockdown_offset + lockdown_transition_period]
result$country_name <- country_adjective

d.ess <- data.frame(as.list(ess(posterior)))
print(d.ess)
colnames(d.ess) <- paste(colnames(d.ess),"ess", sep=".")

d.ess$country = result$country

result <- merge(result, d.ess)

q1 <- data.frame(t(ci(posterior, ci=0.01)),stringsAsFactors=FALSE)

d.medians <- as.list(as.numeric(q1[3,]))
names(d.medians) <- paste(as.vector(q1[1,]), "median", sep=".")
d.medians$country = result$country

result <- merge(result, d.medians)

q95 <- data.frame(t(ci(posterior, ci=0.95)),stringsAsFactors=FALSE)

d.cri95low <- as.list(as.numeric(q95[3,]))
names(d.cri95low) <- paste(as.vector(q95[1,]), "cri95lo", sep=".")
d.cri95low$country = result$country

result <- merge(result, d.cri95low)

d.cri95hi <- as.list(as.numeric(q95[4,]))
names(d.cri95hi) <- paste(as.vector(q95[1,]), "cri95hi", sep=".")
d.cri95hi$country = result$country

result <- merge(result, d.cri95hi)

selection <- 1:dim(posterior)[1]
scount <- length(selection)
draws <- sample(floor(scount/8):scount, quantilePlotSampleSize)
data_sample <- posterior[selection[draws],]

est.died <- data.frame(quantileData(data_sample, function(state, params) { state$died }, 0, lockdown_offset + 200, c(0.05, 0.5, 0.95)))
colnames(est.died) <- c("q5", "q50", "q95")

result$d1.est.died.median <- est.died$q50[lockdown_offset]
result$d1.est.died.cri95lo <- est.died$q5[lockdown_offset]
result$d1.est.died.cri95hi <- est.died$q95[lockdown_offset]

result$d2.est.died.median <- est.died$q50[lockdown_offset + lockdown_transition_period]
result$d2.est.died.cri95lo <- est.died$q5[lockdown_offset + lockdown_transition_period]
result$d2.est.died.cri95hi <- est.died$q95[lockdown_offset + lockdown_transition_period]

result$d2.30.est.died.median <- est.died$q50[lockdown_offset + lockdown_transition_period + 30]
result$d2.30.est.died.cri95lo <- est.died$q5[lockdown_offset + lockdown_transition_period + 30]
result$d2.30.est.died.cri95hi <- est.died$q95[lockdown_offset + lockdown_transition_period + 30]

result$d2.60.est.died.median <- est.died$q50[lockdown_offset + lockdown_transition_period + 60]
result$d2.60.est.died.cri95lo <- est.died$q5[lockdown_offset + lockdown_transition_period + 60]
result$d2.60.est.died.cri95hi <- est.died$q95[lockdown_offset + lockdown_transition_period + 60]

est.Re <- data.frame(quantileData(data_sample, function(state, params) { state$Re }, 0, lockdown_offset + 200, c(0.05, 0.5, 0.95)))
colnames(est.Re) <- c("q5", "q50", "q95")

result$d1.est.Re.median <- est.Re$q50[lockdown_offset]
result$d1.est.Re.cri95lo <- est.Re$q5[lockdown_offset]
result$d1.est.Re.cri95hi <- est.Re$q95[lockdown_offset]

result$d1.30.est.Re.median <- est.Re$q50[lockdown_offset + 30]
result$d1.30.est.Re.cri95lo <- est.Re$q5[lockdown_offset + 30]
result$d1.30.est.Re.cri95hi <- est.Re$q95[lockdown_offset + 30]

result$d1.60.est.Re.median <- est.Re$q50[lockdown_offset + 60]
result$d1.60.est.Re.cri95lo <- est.Re$q5[lockdown_offset + 60]
result$d1.60.est.Re.cri95hi <- est.Re$q95[lockdown_offset + 60]

result$d2.est.Re.median <- est.Re$q50[lockdown_offset + lockdown_transition_period]
result$d2.est.Re.cri95lo <- est.Re$q5[lockdown_offset + lockdown_transition_period]
result$d2.est.Re.cri95hi <- est.Re$q95[lockdown_offset + lockdown_transition_period]

result$d2.30.est.Re.median <- est.Re$q50[lockdown_offset + lockdown_transition_period + 30]
result$d2.30.est.Re.cri95lo <- est.Re$q5[lockdown_offset + lockdown_transition_period + 30]
result$d2.30.est.Re.cri95hi <- est.Re$q95[lockdown_offset + lockdown_transition_period + 30]

result$d2.60.est.Re.median <- est.Re$q50[lockdown_offset + lockdown_transition_period + 60]
result$d2.60.est.Re.cri95lo <- est.Re$q5[lockdown_offset + lockdown_transition_period + 60]
result$d2.60.est.Re.cri95hi <- est.Re$q95[lockdown_offset + lockdown_transition_period + 60]

est.Rt <- data.frame(quantileData(data_sample, function(state, params) { state$Rt }, 10, lockdown_offset + 200, c(0.05, 0.5, 0.95)))
colnames(est.Rt) <- c("q5", "q50", "q95")

result$d1.est.Rt.median <- est.Rt$q50[10 + lockdown_offset]
result$d1.est.Rt.cri95lo <- est.Rt$q5[10 + lockdown_offset]
result$d1.est.Rt.cri95hi <- est.Rt$q95[10 + lockdown_offset]

result$d1.30.est.Rt.median <- est.Rt$q50[10 + lockdown_offset + 30]
result$d1.30.est.Rt.cri95lo <- est.Rt$q5[10 + lockdown_offset + 30]
result$d1.30.est.Rt.cri95hi <- est.Rt$q95[10 + lockdown_offset + 30]

result$d1.60.est.Rt.median <- est.Rt$q50[10 + lockdown_offset + 60]
result$d1.60.est.Rt.cri95lo <- est.Rt$q5[10 + lockdown_offset + 60]
result$d1.60.est.Rt.cri95hi <- est.Rt$q95[10 + lockdown_offset + 60]

result$d2.est.Rt.median <- est.Rt$q50[10 + lockdown_offset + lockdown_transition_period]
result$d2.est.Rt.cri95lo <- est.Rt$q5[10 + lockdown_offset + lockdown_transition_period]
result$d2.est.Rt.cri95hi <- est.Rt$q95[10 + lockdown_offset + lockdown_transition_period]

result$d2.30.est.Rt.median <- est.Rt$q50[10 + lockdown_offset + lockdown_transition_period + 30]
result$d2.30.est.Rt.cri95lo <- est.Rt$q5[10 + lockdown_offset + lockdown_transition_period + 30]
result$d2.30.est.Rt.cri95hi <- est.Rt$q95[10 + lockdown_offset + lockdown_transition_period + 30]

result$d2.60.est.Rt.median <- est.Rt$q50[10 + lockdown_offset + lockdown_transition_period + 60]
result$d2.60.est.Rt.cri95lo <- est.Rt$q5[10 + lockdown_offset + lockdown_transition_period + 60]
result$d2.60.est.Rt.cri95hi <- est.Rt$q95[10 + lockdown_offset + lockdown_transition_period + 60]

result$d1min1.est.Rt.median <- est.Rt$q50[10 + lockdown_offset - 1]
result$d1min1.est.Rt.cri95lo <- est.Rt$q5[10 + lockdown_offset - 1]
result$d1min1.est.Rt.cri95hi <- est.Rt$q95[10 + lockdown_offset - 1]

result$d1min2.est.Rt.median <- est.Rt$q50[10 + lockdown_offset - 2]
result$d1min2.est.Rt.cri95lo <- est.Rt$q5[10 + lockdown_offset - 2]
result$d1min2.est.Rt.cri95hi <- est.Rt$q95[10 + lockdown_offset - 2]

result$d1min3.est.Rt.median <- est.Rt$q50[10 + lockdown_offset - 3]
result$d1min3.est.Rt.cri95lo <- est.Rt$q5[10 + lockdown_offset - 3]
result$d1min3.est.Rt.cri95hi <- est.Rt$q95[10 + lockdown_offset - 3]

est.fr0Rt <- data.frame(quantileData(data_sample, function(state, params) { state$Rt / params[1] }, 10, lockdown_offset + 200, c(0.05, 0.5, 0.95)))
colnames(est.fr0Rt) <- c("q5", "q50", "q95")

result$d1.est.frRt0.median <- est.fr0Rt$q50[10 + lockdown_offset]
result$d1.est.frRt0.cri95lo <- est.fr0Rt$q5[10 + lockdown_offset]
result$d1.est.frRt0.cri95hi <- est.fr0Rt$q95[10 + lockdown_offset]

result$d1min1.est.frRt0.median <- est.fr0Rt$q50[10 + lockdown_offset - 1]
result$d1min1.est.frRt0.cri95lo <- est.fr0Rt$q5[10 + lockdown_offset - 1]
result$d1min1.est.frRt0.cri95hi <- est.fr0Rt$q95[10 + lockdown_offset - 1]

result$d1min2.est.frRt0.median <- est.fr0Rt$q50[10 + lockdown_offset - 2]
result$d1min2.est.frRt0.cri95lo <- est.fr0Rt$q5[10 + lockdown_offset - 2]
result$d1min2.est.frRt0.cri95hi <- est.fr0Rt$q95[10 + lockdown_offset - 2]

result$d1min3.est.frRt0.median <- est.fr0Rt$q50[10 + lockdown_offset - 3]
result$d1min3.est.frRt0.cri95lo <- est.fr0Rt$q5[10 + lockdown_offset - 3]
result$d1min3.est.frRt0.cri95hi <- est.fr0Rt$q95[10 + lockdown_offset - 3]

## Rt2 / Rt@d1
est.fr1Rt <- data.frame(quantileData(data_sample, function(state, params) { params[3] / state$Rt }, 10, lockdown_offset + 200, c(0.05, 0.5, 0.95)))
colnames(est.fr1Rt) <- c("q5", "q50", "q95")

result$d1.est.frRt2.median <- est.fr1Rt$q50[10 + lockdown_offset]
result$d1.est.frRt2.cri95lo <- est.fr1Rt$q5[10 + lockdown_offset]
result$d1.est.frRt2.cri95hi <- est.fr1Rt$q95[10 + lockdown_offset]

result$d1min1.est.frRt2.median <- est.fr1Rt$q50[10 + lockdown_offset - 1]
result$d1min1.est.frRt2.cri95lo <- est.fr1Rt$q5[10 + lockdown_offset - 1]
result$d1min1.est.frRt2.cri95hi <- est.fr1Rt$q95[10 + lockdown_offset - 1]

result$d1min2.est.frRt2.median <- est.fr1Rt$q50[10 + lockdown_offset - 2]
result$d1min2.est.frRt2.cri95lo <- est.fr1Rt$q5[10 + lockdown_offset - 2]
result$d1min2.est.frRt2.cri95hi <- est.fr1Rt$q95[10 + lockdown_offset - 2]

result$d1min3.est.frRt2.median <- est.fr1Rt$q50[10 + lockdown_offset - 3]
result$d1min3.est.frRt2.cri95lo <- est.fr1Rt$q5[10 + lockdown_offset - 3]
result$d1min3.est.frRt2.cri95hi <- est.fr1Rt$q95[10 + lockdown_offset - 3]


est.contribLD <- data.frame(quantileData(data_sample, function(state, params) { (state$Rt - params[3]) / (params[1] - params[3]) }, 10, lockdown_offset + 200, c(0.05, 0.5, 0.95)))
colnames(est.contribLD) <- c("q5", "q50", "q95")

result$est.contribLD.median <- est.contribLD$q50[10 + lockdown_offset]
result$est.contribLD.cri95lo <- est.contribLD$q5[10 + lockdown_offset]
result$est.contribLD.cri95hi <- est.contribLD$q95[10 + lockdown_offset]

result$est.min1.contribLD.median <- est.contribLD$q50[10 + lockdown_offset - 1]
result$est.min1.contribLD.cri95lo <- est.contribLD$q5[10 + lockdown_offset - 1]
result$est.min1.contribLD.cri95hi <- est.contribLD$q95[10 + lockdown_offset - 1]

result$est.min2.contribLD.median <- est.contribLD$q50[10 + lockdown_offset - 2]
result$est.min2.contribLD.cri95lo <- est.contribLD$q5[10 + lockdown_offset - 2]
result$est.min2.contribLD.cri95hi <- est.contribLD$q95[10 + lockdown_offset - 2]

result$est.min3.contribLD.median <- est.contribLD$q50[10 + lockdown_offset - 3]
result$est.min3.contribLD.cri95lo <- est.contribLD$q5[10 + lockdown_offset - 3]
result$est.min3.contribLD.cri95hi <- est.contribLD$q95[10 + lockdown_offset - 3]

print(c(result$d1.est.frRt2.median, result$d1.est.frRt0.median, result$est.contribLD.median))

## offset

integrate <- function(sample, fun)
{
    sampleSize <- dim(sample)[1]

    result <- c()

    for (i in 1:sampleSize) {
        params <- sample[i,]
        params <- transformParams(unlist(params, use.names=FALSE))
        state <- calculateModel(params, lockdown_offset + 100)
        result <- c(result, fun(state))
    }

    result
}

offsets <- quantile(integrate(data_sample, function(state, params) { state$offset }), c(0.05, 0.5, 0.95, 0))

result$p0.d1.median <- offsets[2] + lockdown_offset
result$p0.d1.cri95lo <- offsets[1] + lockdown_offset
result$p0.d1.cri95hi <- offsets[3] + lockdown_offset

quantileAllData <- function(sample, fun, period, quantiles, offset)
{
    simperiod = period * 3
    
    sampleSize <- dim(sample)[1]
    data <- matrix(nrow=(simperiod + offset), ncol=sampleSize)

    for (i in 1:sampleSize) {
        params <- sample[i,]
        params <- transformParams(unlist(params, use.names=FALSE))
        state <- calculateModel(params, simperiod)
        v <- takeAndPad(fun(state), state$offset - offset, simperiod + offset)
        data[,i] = v
    }

    result <- data.frame()
    result <- matrix(nrow=period, ncol=length(quantiles))
    for (j in 1:(dim(result)[1])) {
        result[j,] = quantile(data[j,], quantiles, na.rm=T)
    }
    result
}

max_offset <- offsets[4]

infected <- data.frame(quantileAllData(data_sample, function(state, params) {
    if ("In" %in% names(state)) {
        return(state$E + state$In + state$Is)
    } else {
        return(state$E + state$I)
    } }, lockdown_offset + 200, c(0.05, 0.5, 0.95), max_offset))
colnames(infected) <- c("q5", "q50", "q95")

result$d1.infected.median <- infected$q50[max_offset + lockdown_offset]
result$d1.infected.cri95lo <- infected$q5[max_offset + lockdown_offset]
result$d1.infected.cri95hi <- infected$q95[max_offset + lockdown_offset]

result$d2.infected.median <- infected$q50[max_offset + lockdown_offset + lockdown_transition_period]
result$d2.infected.cri95lo <- infected$q5[max_offset + lockdown_offset + lockdown_transition_period]
result$d2.infected.cri95hi <- infected$q95[max_offset + lockdown_offset + lockdown_transition_period]

if (max_offset + lockdown_offset - 10 > 0) {
    result$d1.min10.infected.median <- infected$q50[max_offset + lockdown_offset - 10]
    result$d1.min10.infected.cri95lo <- infected$q5[max_offset + lockdown_offset - 10]
    result$d1.min10.infected.cri95hi <- infected$q95[max_offset + lockdown_offset - 10]
} else {
    result$d1.min10.infected.median <- 0
    result$d1.min10.infected.cri95lo <- 0
    result$d1.min10.infected.cri95hi <- 0
}

if (max_offset + lockdown_offset - 20 > 0) {
    result$d1.min20.infected.median <- infected$q50[max_offset + lockdown_offset - 20]
    result$d1.min20.infected.cri95lo <- infected$q5[max_offset + lockdown_offset - 20]
    result$d1.min20.infected.cri95hi <- infected$q95[max_offset + lockdown_offset - 20]
} else {
    result$d1.min20.infected.median <- 0
    result$d1.min20.infected.cri95lo <- 0
    result$d1.min20.infected.cri95hi <- 0
}

if (max_offset + lockdown_offset - 30 > 0) {
    result$d1.min30.infected.median <- infected$q50[max_offset + lockdown_offset - 30]
    result$d1.min30.infected.cri95lo <- infected$q5[max_offset + lockdown_offset - 30]
    result$d1.min30.infected.cri95hi <- infected$q95[max_offset + lockdown_offset - 30]
} else {
    result$d1.min30.infected.median <- 0
    result$d1.min30.infected.cri95lo <- 0
    result$d1.min30.infected.cri95hi <- 0
}

if (max_offset + lockdown_offset - 60 > 0) {
    result$d1.min60.infected.median <- infected$q50[max_offset + lockdown_offset - 60]
    result$d1.min60.infected.cri95lo <- infected$q5[max_offset + lockdown_offset - 60]
    result$d1.min60.infected.cri95hi <- infected$q95[max_offset + lockdown_offset - 60]
} else {
    result$d1.min60.infected.median <- 0
    result$d1.min60.infected.cri95lo <- 0
    result$d1.min60.infected.cri95hi <- 0
}

if (max_offset + lockdown_offset - 90 > 0) {
    result$d1.min90.infected.median <- infected$q50[max_offset + lockdown_offset - 90]
    result$d1.min90.infected.cri95lo <- infected$q5[max_offset + lockdown_offset - 90]
    result$d1.min90.infected.cri95hi <- infected$q95[max_offset + lockdown_offset - 90]
} else {
    result$d1.min90.infected.median <- 0
    result$d1.min90.infected.cri95lo <- 0
    result$d1.min90.infected.cri95hi <- 0
}

## Analyze predictive performance
est.deadi <- data.frame(quantileData(data_sample, function(state, params) { state$deadi }, 0, lockdown_offset + 200, c(0.05, 0.5, 0.95)))
colnames(est.deadi) <- c("q5", "q50", "q95")

fitl <-length(dmorti) - evaluation_data_count

result$fit.logl <- sum(dnbinom(dmorti[1:fitl],
                               mu=pmax(1E-10, est.deadi$q50[1:fitl]), size=mort_nbinom_size, log=T))

if (evaluation_data_count > 0) {
    result$eval.logl <- sum(dnbinom(dmorti[(fitl + 1):length(dmorti)],
                                    mu=pmax(1E-10, est.deadi$q50[(fitl + 1):length(dmorti)]),
                                    size=mort_nbinom_size, log=T))
} else {
    result$eval.logl <- 0
}

result$fit.rmse <- sqrt(sum((dmorti[1:fitl] - est.deadi$q50[1:fitl])^2))

if (evaluation_data_count > 0) {
    result$eval.rmse <- sqrt(sum((dmorti[(fitl + 1):length(dmorti)] - est.deadi$q50[(fitl + 1):length(dmorti)])^2))
} else {
    result$eval.rmse <- 0
}

gm <- c

gmbl <- subset(gm[1:o$par[1],])
gmal <- subset(gm[o$par[2]:(dim(gm)[1]),])

result$gm.bl.retail_and_recreation_percent_change_from_baseline <- mean(gmbl$retail_and_recreation_percent_change_from_baseline, na.rm = T)
result$gm.al.retail_and_recreation_percent_change_from_baseline <- mean(gmal$retail_and_recreation_percent_change_from_baseline, na.rm = T)

result$gm.bl.grocery_and_pharmacy_percent_change_from_baseline <- mean(gmbl$grocery_and_pharmacy_percent_change_from_baseline, na.rm = T)
result$gm.al.grocery_and_pharmacy_percent_change_from_baseline <- mean(gmal$grocery_and_pharmacy_percent_change_from_baseline, na.rm = T)

result$gm.bl.parks_percent_change_from_baseline <- mean(gmbl$parks_percent_change_from_baseline, na.rm = T)
result$gm.al.parks_percent_change_from_baseline <- mean(gmal$parks_percent_change_from_baseline, na.rm = T)

result$gm.bl.transit_stations_percent_change_from_baseline <- mean(gmbl$transit_stations_percent_change_from_baseline, na.rm = T)
result$gm.al.transit_stations_percent_change_from_baseline <- mean(gmal$transit_stations_percent_change_from_baseline, na.rm = T)

result$gm.bl.workplaces_percent_change_from_baseline <- mean(gmbl$workplaces_percent_change_from_baseline, na.rm = T)
result$gm.al.workplaces_percent_change_from_baseline <- mean(gmal$workplaces_percent_change_from_baseline, na.rm = T)

result$gm.bl.residential_percent_change_from_baseline <- mean(gmbl$residential_percent_change_from_baseline, na.rm = T)
result$gm.al.residential_percent_change_from_baseline <- mean(gmal$residential_percent_change_from_baseline, na.rm = T)

## Must be last since it modifies the data
dic = DIC(posterior1)

print(dic)

result$DIC = dic$DIC
result$DIC2 = dic$DIC2
result$IC = dic$IC
result$pD = dic$pD
result$pV = dic$pV
result$Dbar = dic$Dbar
result$Dhat = dic$Dhat

write.csv(result, paste(country2,"analysis.csv",sep='_'))

