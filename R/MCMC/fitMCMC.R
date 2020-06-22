require(adaptMCMC)

source("settings.R")

source(data, chdir=T)
source(fitmodel, chdir=T)

source(paste(Rdir, "lib/libfit.R", sep=""))

options(scipen=999)

calcloglMCMC <- function(params) {
    tr <- transformParams(params)
    return(calclogp(tr) + calclogl(tr, NULL))
}

## source it again so that you can override things
source("settings.R")

print(fit.paramnames)

print(c("initial logl: ", calcloglMCMC(init)))

it <- 0
out <<- MCMC(calcloglMCMC, n=1E3, init=init, scale=scales, adapt=200000, acc.rate=0.234)
batch <- out$samples[seq(1, nrow(out$samples), 1E2), ]
colnames(batch) <- fit.paramnames

##
## Do iterations
##
it <- 0

for (i in 1:(72*20)) {
    print(c("Iteration", i))
    out <<- MCMC.add.samples(out, n=1E4)
    batch <- data.frame(out$samples[seq(1, nrow(out$samples), 1E2), ])
    colnames(batch) <- fit.paramnames
    plot(ts(subset(batch, select=fitkeyparamnames)))
    write.csv(batch, file=outputfile)
}
