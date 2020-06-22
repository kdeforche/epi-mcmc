require(drjacoby)

source("settings.R")

source(data, chdir=T)
source(fitmodel, chdir=T)

source(paste(Rdir, "lib/libfit.R", sep=""))

options(scipen=999)

calcloglMCMC <- function(params, x) {
    return(calclogl(transformParams(params)))
}

## source it again so that you can override things
source("settings.R")

print(fit.paramnames)

print(c("initial logl: ", calcloglMCMC(init), " logl prior: ", calclogp(transformParams(init))))

for (chain in 1:4) {
  r_mcmc_out <- run_mcmc(data = dhospi,
                       df_params = df_params,
                       loglike = calcloglMCMC,
                       logprior = calclogp,
                       burnin = 2e3,
                       samples = 2e3,
                       rungs = 20,
                       chains = 1,
                       GTI_pow = 1)
  print(r_mcmc_out$diagnostics$ess)
  print(data.frame(r_mcmc_out$diagnostics$mc_accept))
  r_mcmc_out$output$chain <- paste("chain",chain,sep='')
  chain_outputfile <- paste(outputfile,"_",chain,sep='')
  write.csv(subset(r_mcmc_out$output, rung=="rung1" & stage=="sampling"), file=chain_outputfile)
  write.csv(r_mcmc_out$output, file=paste("all_", chain_outputfile, sep=''))
}

## r_mcmc_out <- run_mcmc(data = dhospi,
##                        df_params = df_params,
##                        loglike = calcloglMCMC,
##                        logprior = calclogp,
##                        burnin = 1e3,
##                        samples = 1e3,
##                        rungs = 20,
##                        chains = 1,
##                        GTI_pow = 1)

## print(r_mcmc_out$diagnostics$ess)
## print(data.frame(r_mcmc_out$diagnostics$mc_accept))
## write.csv(subset(r_mcmc_out$output, rung=="rung1" & stage=="sampling"), file=outputfile)
## write.csv(r_mcmc_out$output, file=paste("all_", outputfile, sep=''))
