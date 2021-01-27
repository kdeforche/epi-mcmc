require(drjacoby)

source("settings.R")

source(data, chdir=T)
source(fitmodel, chdir=T)

if (file.exists("max.csv")) {
    print("Loading initial state from max.csv")
    r <- read.csv("max.csv")
    init <- r$x[2:(length(init) + 1)]
    df_params$init = init
}

source(paste(Rdir, "lib/libfit.R", sep=""))

options(scipen=999)

calcloglMCMC <- function(params, x, misc) {
    return(calclogl(transformParams(params)))
}

## source it again so that you can override things
source("settings.R")

print(fit.paramnames)

print(c("initial logl: ", calcloglMCMC(init), " logl prior: ", calclogp(transformParams(init))))

cl <- NULL
chains <- 1
## cores <- parallel::detectCores()
## cl <- parallel::makeCluster(cores)

## parallel::clusterCall(cl, function() {
##     source("settings.R")
##     source(data, chdir=T)
##     source(fitmodel, chdir=T)
##     source(paste(Rdir, "lib/libfit.R", sep=""))
##     it <<- 0
## })

data = list(x=2)

for (chain in 1:4) {
  r_mcmc_out <- run_mcmc(data = data,
                       df_params = df_params,
                       loglike = calcloglMCMC,
                       logprior = calclogp,
                       burnin = 1e3,
                       samples = 0.5e3,
                       rungs = 20,
                       chains = chains,
                       cluster = cl,
                       GTI_pow = 1)
  print(r_mcmc_out$diagnostics$ess)
  print(data.frame(r_mcmc_out$diagnostics$mc_accept))
  r_mcmc_out$output$chain <- paste("chain",chain,sep='')
  chain_outputfile <- paste(outputfile,"_",chain,sep='')
  write.csv(subset(r_mcmc_out$output, rung==20 & phase=="sampling"), file=chain_outputfile)
##  write.csv(r_mcmc_out$output, file=paste("all_", chain_outputfile, sep=''))
}
