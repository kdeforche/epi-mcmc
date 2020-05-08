library("bayesplot")
library("ggplot2")
require(mcmcse)
require(bayestestR)

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

system(paste("mkdir ", outputdir))

read <- function() {
    all <- readData(inputfiles)

    posterior <- all

    print(ess(posterior))

    plot(ts(subset(posterior, select=keyparamnames)))

    ## compute credibility intervals
    print(data.frame(ci(posterior, ci=0.01)))
    print(data.frame(ci(posterior, ci=0.50)))
    print(data.frame(ci(posterior, ci=0.95)))

    posterior
}

densityPlot <- function() {
    ## sample a subset of data to show in density plots
    selection <- which(posterior[,"Tinc"] < 100)
    scount <- length(selection)
    draws <- sample(floor(scount/8):scount, densityPlotSampleSize)
    display_sample <- posterior[selection[draws],][0:-1]

    ## show daily incidence plots + posterior density
    predict_daily_plot(display_sample, F, NULL)

    ## alternative, only hospital incidence
    ##predict_daily_plot(display_sample, F, alpha("red", 0.1))
    ##title('New hospitalisations per outputdir')
    ##legend("topleft", inset=0.02, legend=c("Hospitalisations"),
    ##       col=c("red"),lty=1)

    ##predict_daily_plot(display_sample, F, alpha("#CCCCCC", 0.1))
    ##predict_daily_plot(display_sample, T, alpha("#000066", 0.3))
    ##predict_daily_plot(display_sample, T, alpha("#660000", 0.3))

    #png(paste(outputdir, "/forecast-time.png", sep=""), width=500, height=500)
    #predict_daily_plot(display_sample, F, NULL)
    #dev.off()

    #pdf(paste(outputdir, "/forecast-time.pdf", sep=""), width=8, height=8)
    #predict_daily_plot(display_sample, F, NULL)
    #dev.off()
}

calcloglMCMC <- function(params) {
    params <- transformParams(unlist(params, use.names=FALSE))
    return(calclogl(params))
}

posterior <- read()
densityPlot()

posterior$logl <- numeric(dim(posterior)[1])
for (i in 1:dim(posterior)[1]) {
    params <- posterior[i,0:-1]
    posterior$logl[i] <- calcloglMCMC(params)
}

#
# from here on only totally random snippets
#
mcmc_hist(posterior, pars = keyparamnames)

#mcmc_hist(posterior, pars = c("IFR"), binwidth=0.0001)
#mcmc_hist(posterior, pars = c("Beta"))
#mcmc_hist(posterior, pars = c("HL", "DL"), binwidth=0.1)
#mcmc_hist(posterior, pars = c("WZC"))

png("Ne.png", width=500, height=500)
mcmc_hist(posterior, pars = c("Ne"))
dev.off()

png("Rt.png", width=500, height=500)
pdf(paste(outputdir, "/Rt.pdf", sep=""), width=6, height=6)
mcmc_intervals(posterior, pars=rev(c("y.R0", "o.R0", "yo.R0", "y.Rt", "o.Rt", "yo.Rt")),
               prob = 0.5, prob_outer = 0.95) +
    scale_x_continuous(name="Reproduction number", breaks=1:6) +
    scale_y_discrete("", labels=rev(c("R0,y", "R0,o", "R0,yo", "Rt,y", "Rt,o", "Rt,yo"))) +
    geom_vline(xintercept=1.0, linetype="dashed", color='red', size=1.3)
dev.off()

mcmc_areas(posterior, pars = c("R0"), prob = 0.8)

mcmc_areas(posterior,
           pars = c("y.R0", "y.Rt"),
           prob = 0.8)

png("Rt.png", width=500, height=500)
pdf(paste(outputdir, "/Rt.pdf", sep=""), width=4, height=4)
mcmc_areas(posterior,
           pars = c("Rt"),
           prob = 0.8) + scale_x_continuous(name="", limits=c(0, 3))
dev.off()

mcmc_areas(posterior,
           pars = c("Tinc", "Tinf"),
           prob = 0.8) + scale_x_continuous(name="", limits=c(0, 20))

s <- seq(0.00001, 0.035, 0.00001)
ifr_prior1 <- dnorm(s, 0.007, 0.0025)
ifr_prior2 <- dbeta(s, 4.8, 722.69)
ifr_prior3 <- dbeta(s, 2.697931, 406.0796)
ifr_prior4 <- dbeta(s, 1.7243, 259.53)
ifr_prior5 <- dbeta(s, 10.8, 1627)
ifr_prior6 <- dbeta(s, 2.24, 744)
ifr_prior7 <- dbeta(s, 0.994, 330)
y.ifr_prior <- dbeta(s, 1.263586, 1402.721)
o.ifr_prior <- dbeta(s, 1.31346, 54.81763)
plot(s, y.ifr_prior, type='l')
lines(s, o.ifr_prior, type='l')
lines(s, ifr_prior6, type='l', col='red')
lines(s, ifr_prior5, type='l', col='blue')
lines(s, ifr_prior4, type='l', col='green')

mcmc_areas(posterior,
           pars = c("y.IFR"),
           prob = 0.8) + scale_x_continuous(name="y.IFR", limits=c(0, 0.015))
#lines(s, ifr_prior1, type='l', col='red')
#lines(s, ifr_prior2, type='l', col='blue')
#lines(s, ifr_prior3, type='l', col='gray')
#lines(s, ifr_prior4, type='l', col='black')
lines(s, ifr_prior5, type='l', col='blue')


## extra plot met:
## - Rt evolutie
## - Rt evolutie
