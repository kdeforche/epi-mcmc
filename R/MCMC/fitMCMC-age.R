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

init <- c(3.6, 0.7, 3.04, 0.12, 0.27, 0.09,
          log(0.01), log(0.3), 14, 13, 1.3, 7, total_deaths_at_lockdown)
scales <- c(0.15, 0.05, 0.15, 0.05, 0.15, 0.05,
            0.05, 0.1, 1, 1, 0.3, 0.3, total_deaths_at_lockdown / 20)
scale <- 0.15 * scales

source("control.R")

print(c("initial logl: ", calcloglMCMC(init)))

##
## Do a few iterations to check that things are setup properly.
##
it <- 0
out <<- metrop(calcloglMCMC, init, scale=scale, nbatch=100, blen=1)
colnames(out$batch) <- cn
out$accept
#graphs()

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
    out <<- metrop(out, init=out$final, scale=scale, nbatch=100, blen=100)
    print(paste(c("============= iteration "), i, " acceptance: ", out$accept))
    colnames(out$batch) <- cn
    ##plot(ts(out$batch))
    write.table(out$batch, file=outputfile, append=T, quote=F, sep=",", col.name=F)

    source("control.R")
}

#iterate()
