require(mcmc)

source("control.R")

source(settings)
source(data, chdir=T)
source(fitmodel, chdir=T)

source(paste(Rdir, "lib/libfit.R", sep=""))

options(scipen=999)

calcloglMCMC <- function(params) {
    return(calclogl(transformParams(params)))
}

cn <- fit.paramnames

init <- c(2.3, 0.5, log(0.02), log(0.3), 10, 9, 2, 5, total_deaths_at_lockdown, 0)
scales <- c(0.15, 0.05, 0.05, 0.1, 1, 1, 0.3, 0.3, total_deaths_at_lockdown / 20, 0.1)
scale <- 0.15 * scales

source("control.R")

print(c("initial logl: ", calcloglMCMC(init)))

##
## Do a few iterations to check that things are setup properly.
##
it <- 0
out <- metrop(calcloglMCMC, init, scale=scale, nbatch=100, blen=1)
colnames(out$batch) <- cn
out$accept
graphs()

##
## Do iterations
##
it <- 0

if (truncate) {
    write.csv(out$batch, file=outputfile)
}

## almost endless loop to calculate
scale <- 0.125 * scales

for (i in 1:10000) {
    out <- metrop(out, init=out$final, scale=scale, nbatch=100, blen=100)
    print(paste(c("============= iteration "), i, " acceptance: ", out$accept))
    colnames(out$batch) <- cn
    plot(ts(out$batch))
    write.table(out$batch, file=outputfile, append=T, quote=F, sep=",", col.name=F)
    source("control.R")
}
