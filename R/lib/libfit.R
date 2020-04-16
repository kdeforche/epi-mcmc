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

##
## Graphs used to monitor the MCMC runs
##
graphs <- function() {
    par(mfrow=c(1,2))

    print(names(state))
    
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

    period <- length(state$hospi)
    len <- period - state$offset + 1
    plot(days[1:len],state$hospi[state$offset:period], type='l', col='red',
         xlab='Date', ylab='Count',
         main='New hospitalisations and deaths per day')

    if ("y.hospi" %in% names(state)) {
        lines(days[1:len],state$y.hospi[state$offset:period], type='l', lty=2, col='red')
        lines(days[1:len],state$o.hospi[state$offset:period], type='l', lty=3, col='red')
        lines(days[1:len],state$y.deadi[state$offset:period], type='l', lty=2, col='blue')
        lines(days[1:len],state$o.deadi[state$offset:period], type='l', lty=3, col='blue')
    } else {
        lines(days[1:len],state$deadi[state$offset:period], type='l', col='blue')
    }

    if (exists("y.dmorti")) {
        points(days[1:length(dmorti)],y.dmorti,col=c("blue"))
        points(days[1:length(dmorti)],o.dmorti,col=c("blue"))
        points(days[1:length(dhospi)],y.dhospi,col=c("red"))
        points(days[1:length(dhospi)],o.dhospi,col=c("red"))
    } else {
        points(days[1:length(dmorti)],dmorti,col=c("blue"))
        points(days[1:length(dhospi)],dhospi,col=c("red"))
    }

    legend("topright", inset=0.02, legend=c("Hospitalisations", "Deaths"),
	   col=c("red", "blue"),lty=1)

    if ("R" %in% names(state)) {
        print(paste("% immune: ", state$R[period]/y.N))
    } else {
        print(paste("% y immune: ", state$y.R[period]/y.N))
        print(paste("% o immune: ", state$o.R[period]/o.N))
    }
    print(paste("deaths: ", state$died[length(state$died) - 2]))
    print(paste("offset: ", state$offset - state$padding))
}
	
