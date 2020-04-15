library("bayesplot")
library("ggplot2")
require(mcmcse)
require(bayestestR)

scales <- 1

source("control.R")
source(settings)

source(data, chdir=T)
source(fitmodel, chdir=T)

source(paste(Rdir, "lib/libMCMC.R", sep=""))

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

# sample a subset of data to show in density plots
selection <- which(posterior[,"IFR"] < 3)
scount <- length(selection)
draws <- sample(floor(scount/8):scount, densityPlotSampleSize)
display_sample <- posterior[selection[draws],][0:-1]

# show daily incidence plots + posterior density
predict_daily_plot(display_sample, F, NULL)

# alternative, only hospital incidence
#predict_daily_plot(display_sample, F, alpha("red", 0.1))
#title('New hospitalisations per outputdir')
#legend("topleft", inset=0.02, legend=c("Hospitalisations"),
#       col=c("red"),lty=1)

#predict_daily_plot(display_sample, F, alpha("#CCCCCC", 0.1))
#predict_daily_plot(display_sample, T, alpha("#000066", 0.3))
#predict_daily_plot(display_sample, T, alpha("#660000", 0.3))

png(paste(outputdir, "/forecast-time.png", sep=""), width=500, height=500)
predict_daily_plot(display_sample, F, NULL)
dev.off()

pdf(paste(outputdir, "/forecast-time.pdf", sep=""), width=8, height=8)
predict_daily_plot(display_sample, F, NULL)
dev.off()

#
# from here on only totally random snippets
#

mcmc_hist(posterior, pars = c("y.R0", "y.Rt", "IFR"))

#mcmc_hist(posterior, pars = c("IFR"), binwidth=0.0001)
#mcmc_hist(posterior, pars = c("Beta"))
#mcmc_hist(posterior, pars = c("HL", "DL"), binwidth=0.1)
#mcmc_hist(posterior, pars = c("WZC"))

mcmc_areas(posterior, pars = c("R0"), prob = 0.8)

mcmc_areas(posterior,
           pars = c("Rt"),
           prob = 0.8) + scale_x_continuous(name="", limits=c(0, 3))

png("Rt.png", width=500, height=500)
pdf(paste(outputdir, "/Rt.pdf", sep=""), width=4, height=4)
mcmc_areas(posterior,
           pars = c("Rt"),
           prob = 0.8) + scale_x_continuous(name="", limits=c(0, 3))
dev.off()

mcmc_areas(posterior,
           pars = c("Tinc", "Tinf"),
           prob = 0.8) + scale_x_continuous(name="", limits=c(0, 20))

s <- seq(0.0001, 0.015, 0.0001)
ifr_prior1 <- dnorm(s, 0.007, 0.0025)
ifr_prior2 <- dbeta(s, 4.8, 722.69)
ifr_prior3 <- dbeta(s, 2.697931, 406.0796)
ifr_prior4 <- dbeta(s, 1.7243, 259.53)
ifr_prior5 <- dbeta(s, 10.8, 1627)
plot(s, ifr_prior5, type='l')

mcmc_areas(posterior,
           pars = c("IFR"),
           prob = 0.8) + scale_x_continuous(name="IFR", limits=c(0, 0.015))
#lines(s, ifr_prior1, type='l', col='red')
#lines(s, ifr_prior2, type='l', col='blue')
#lines(s, ifr_prior3, type='l', col='gray')
#lines(s, ifr_prior4, type='l', col='black')
lines(s, ifr_prior5, type='l', col='blue')
