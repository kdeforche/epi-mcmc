require(adaptMCMC)

source("control.R")

source(settings)
source(data, chdir=T)
source(fitmodel, chdir=T)

source(paste(Rdir, "lib/libfit.R", sep=""))

options(scipen=999)

source("control.R")

calcloglMCMC <- function(params) {
    return(calclogl(transformParams(params)))
}

cn <- fit.paramnames

init <- c(3.6, 0.7, 3.04, 0.12, 0.27, 0.09,
          log(0.01), log(0.3), log(1),
          14, 13, 14, 13, 1.3, 7, total_deaths_at_lockdown, log(25))
scales <- c(0.15, 0.05, 0.15, 0.05, 0.15, 0.05,
            0.05, 0.05, 0.05,
            1, 1, 1, 1, 0.3, 0.3, total_deaths_at_lockdown / 20, 0.1)

print(c("initial logl: ", calcloglMCMC(init)))

it <- 0
out <<- MCMC(calcloglMCMC, n=1E3, init=init, scale=scales, adapt=T, acc.rate=0.234)
batch <- out$samples[seq(1, nrow(out$samples), 1E2), ]
colnames(batch) <- cn

##
## Do iterations
##
it <- 0

for (i in 1:10000) {
    out <<- MCMC.add.samples(out, n=1E4)
    batch <- data.frame(out$samples[seq(1, nrow(out$samples), 1E2), ])
    colnames(batch) <- cn
    plot(ts(subset(batch, select=fitkeyparamnames)))
    write.csv(batch, file=outputfile)
}
