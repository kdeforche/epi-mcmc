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
maxLogL <- -1E10

graphs <- function(logl) {
    ## only if a device is open already
    if (.Device == "null device") {
        return(0)
    }
   
    par(mfrow=c(1,2))

    days <- seq(dstartdate, dstartdate + length(dmort) + 30, 1)

    period <- length(state$y.deadi)
    len <- period - state$offset + 1

    if (state$offset < 1 | len < 2) {
        plot(days[1:5], state$y.deadi[1:5], type='l', col='red',
             xlab='Date', ylab='Count',
             main='Cumulative deaths')
        plot(days[1:5], state$y.deadi[1:5], type='l', col='red',
             xlab='Date', ylab='Count',
             main='New cases and deaths per day')
        return (0)
    }

    plot(days[1:len], state$died[state$offset:(state$offset + len - 1)], type='l', col='blue',
         xlab='Date', ylab='Cumulative count',
         main='Cumulative cases and deaths',
         ylim=c(10, 1E6),
         log="y")
    points(days[1:length(dmort)],dmort,col='blue')

    plot(days[1:len],state$y.deadi[state$offset:period], type='l', lty=2, col='red',
         xlab='Date', ylab='Count', ylim=c(0.1, 20000),
         main=paste('Fit, logl=', logl), log="y")

    lines(days[1:len],state$y.casei[state$offset:period], type='l', lty=2, col='darkgreen')
    lines(days[1:len],state$o.casei[state$offset:period], type='l', lty=3, col='darkgreen')
    lines(days[1:len],state$o.deadi[state$offset:period], type='l', lty=3, col='blue')
    lines(days[1:len],(state$y.hospi + state$o.hospi)[state$offset:period], type='l', col='orange')
    lines(days[1:len],(state$y.hospi)[state$offset:period], type='l', lty=2, col='orange')
    lines(days[1:len],(state$o.hospi)[state$offset:period], type='l', lty=3, col='orange')

    points(days[1:length(y.dcasei)],y.dcasei,col=c("darkgreen"))
    points(days[1:length(o.dcasei)],o.dcasei,col=c("darkgreen"))
    points(days[1:length(y.dmorti)],y.dmorti,col=c("red"))
    points(days[1:length(o.dmorti)],o.dmorti,col=c("blue"))
    points(days[1:length(dhospi)],dhospi,col=c("orange"))

    lines(days[1:len],state$Re[state$offset:period], type='l', lty=1, lwd=2, col='black')
}
