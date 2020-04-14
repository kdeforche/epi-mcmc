###############
## Default values for fitting data. You may want to override that in your script.
###############

##
## How to estimate dispersion parameter (size):
## Fit a good model (e.g. from an MCMC batch with low acceptance rate)
##  m = mu, s = data
##  find size so that mean((m - s)^2 / (m + m^2/size)) = 1
##  take size a bit smaller to be safe (because ML fit minimizes variance)
##
## smaller size ("dispersion parameter"), larger variation
## r -> inf : poisson
## r -> 1 : var ~ mu^2
##
hosp_nbinom_size = 30
mort_nbinom_size = 90

##
## Length of model calculations needed to fit the data (you need to take ample time before first
## cases to allow the epidemic to grow.
##
FitTotalPeriod = 80

####################
## Fit library functions
####################

estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

## Global variable that stores calculated model data
state <- NULL

## Iteration counter to log periodically during log likelihood calculations
it <- 0

calclogl <- function(params) {
    beta0 <- params[1]
    betat <- params[2]
    hosp_rate <- params[3]
    died_rate <- params[4]
    hosp_latency <- params[5]
    died_latency <- params[6]
    Tinf <- params[7]
    Tinc <- params[8]
    mort_lockdown_threshold <- params[9]
    hosp_rate_change <- params[10]

    if (beta0 < 0.01) {
        print(paste("invalid beta0", beta0))
        return(-Inf)
    }

    if (betat < 0.01) {
        print(paste("invalid betat", betat))
        return(-Inf)
    }

    if (hosp_rate < 0.0001) {
        print(paste("invalid hosp_rate", hosp_rate))
        return(-Inf)
    }

    if (died_rate < 0.0001 && betat > 0.2) {
        print(paste("invalid died_rate", died_rate))
        return(-Inf)
    }

    if (hosp_latency < 2 || hosp_latency > 30) {
        print(paste("invalid hosp_latency", hosp_latency))
        return(-Inf)
    }

    if (died_latency < 2 || died_latency > 30) {
        ## print(paste("invalid died_latency", died_latency))
        return(-Inf)
    }

    if (Tinf < 1.1 || Tinf > 20) {
        print(paste("invalid Tinf", Tinf))
        return(-Inf)
    }

    if (Tinc < 1.1 || Tinc > 20) {
        print(paste("invalid Tinc", Tinc))
        return(-Inf)
    }

    if (hosp_rate_change < 1.0) {
        print(paste("invalid hosp_rate_change", hosp_rate_change))
        return(-Inf)
    }

    state <<- calculateModel(params, FitTotalPeriod)

    if (state$offset == InvalidDataOffset)
       return(-Inf)

    logl <- 0

    #R0 <- beta0 * Tinf

    #logl <- logl + dnorm(R0, mean=6, sd=0.5, log=T)
    
    #logl <- logl + dnorm(beta0, mean=0.1, sd=3, log=T)
    #logl <- logl + dnorm(betat, mean=0.1, sd=3, log=T)
    #logl <- logl + dgamma(beta0, 1, 1, log=T)
    #logl <- logl + dgamma(betat, 1, 1, log=T)
    #logl <- logl + dgamma(hosp_rate, shape=0.01, rate=1, log=T)
    
    #logl <- logl + dbeta(died_rate, 10.8, 1627, log=T) # estBetaParams(0.0066, 0.002^2)
    #logl <- logl + dbeta(died_rate, 4.8, 722.69, log=T) # estBetaParams(0.0066, 0.003^2)
    #logl <- logl + dbeta(died_rate, 2.697931, 406.0796, log=T) # estBetaParams(0.0066, 0.004^2)
    logl <- logl + dbeta(died_rate, 1.7243, 259.53, log=T) # estBetaParams(0.0066, 0.005^2)
    #logl <- logl + dbeta(died_rate, 0.9, 0.9, log=T) # https://stats.stackexchange.com/questions/297901/choosing-between-uninformative-beta-priors
    logl <- logl + dnorm(hosp_latency, mean=10, sd=20, log=T)
    logl <- logl + dnorm(died_latency, mean=10, sd=20, log=T)
    logl <- logl + dnorm(Tinc, mean=5, sd=3, log=T)
    logl <- logl + dnorm(Tinf, mean=5, sd=3, log=T)
    logl <- logl + dnbinom(total_deaths_at_lockdown, mu=pmax(0.1, mort_lockdown_threshold), size=mort_nbinom_size, log=T)
    logl <- logl + dnorm(hosp_rate_change, mean=1, sd=0.3, log=T)

    dstart <- state$offset
    dend <- state$offset + length(dhospi) - 1

    if (dend > length(state$hospi)) {
       print("=========================== Increase FitTotalPeriod ===================")
       return(-Inf)
    }

    dh <- dhospi
    loglA <- sum(dnbinom(dh, mu=pmax(0.1, state$hospi[dstart:dend]), size=hosp_nbinom_size, log=T))

    dend <- state$offset + length(dmorti) - 1

    if (dend > length(state$deadi)) {
       print("=========================== Increase FitTotalPeriod ===================")
       return(-Inf)
    }

    loglB <- sum(dnbinom(dmorti, mu=pmax(0.1, state$deadi[dstart:dend] / death_underreporting_factor), size=mort_nbinom_size, log=T))

    it <<- it + 1
    if (it %% 1000 == 0) {
        print(params)
	print(c(it, logl + loglA + loglB))
	graphs()
    }

    logl + loglA + loglB
}

##
## Graphs used to monitor the MCMC runs
##
graphs <- function() {
    par(mfrow=c(1,2))

    days <- seq(dstartdate, dstartdate + length(dhosp) + 30, 1)

    len <- 40
    plot(days[1:len], state$hosp[state$offset:(state$offset + len - 1)], type='l', col='red',
         ylim=c(0, 10000),
         xlab='Date', ylab='Cumulative count',
         main='Cumulative hospitalisations and deaths')
    lines(days[1:len],state$died[state$offset:(state$offset + len - 1)], type='l', col='blue')
    points(days[1:length(dhospi)],dhosp,col='red')
    points(days[1:length(dmort)],dmort,col='blue')
    legend("topleft", inset=0.02, legend=c("Hospitalisations", "Deaths"),
	   col=c("red", "blue"),lty=1)

    period <- length(state$R)
    len <- period - state$offset + 1
    plot(days[1:len],state$hospi[state$offset:period], type='l', col='red',
         xlab='Date', ylab='Count',
         main='New hospitalisations and deaths per day')
    lines(days[1:len],state$deadi[state$offset:period], type='l', col='blue')
    points(days[1:length(dhospi)],dhospi,col=c("red"))
    points(days[1:length(dmorti)],dmorti,col=c("blue"))
    legend("topright", inset=0.02, legend=c("Hospitalisations", "Deaths"),
	   col=c("red", "blue"),lty=1)

    print(paste("% immune: ", state$R[period]/N))
    print(paste("deaths: ", state$died[length(state$died) - 2]))
    print(paste("offset: ", state$offset - state$padding))
}
	
