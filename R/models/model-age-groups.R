source("model.R")
require(Rcpp)

cppFunction('NumericVector odesimstepc(double yS, double yE, double yI, double yR, double oS, double oE, double oI, double oR, double ybeta, double obeta, double yobeta, double a, double gamma, int yN, int oN) {

   double ygot_infected = std::min(1.0, ybeta * yI / yN) * yS;
   double ygot_infectious = a * yE;
   double ygot_removed = gamma * yI;

   double ydeltaS = -ygot_infected;
   double ydeltaE = ygot_infected - ygot_infectious;
   double ydeltaI = ygot_infectious - ygot_removed;
   double ydeltaR = ygot_removed;

   double ogot_infected = std::min(1.0, obeta * oI / oN + yobeta * yI / yN) * oS;
   double ogot_infectious = a * oE;
   double ogot_removed = gamma * oI;

   double odeltaS = -ogot_infected;
   double odeltaE = ogot_infected - ogot_infectious;
   double odeltaI = ogot_infectious - ogot_removed;
   double odeltaR = ogot_removed;

   NumericVector out(8);

   out[0] = ydeltaS;
   out[1] = ydeltaE;
   out[2] = ydeltaI;
   out[3] = ydeltaR;

   out[4] = odeltaS;
   out[5] = odeltaE;
   out[6] = odeltaI;
   out[7] = odeltaR;

   return out;
}')

ode.simstep <- function(i, state.seir, parameters) {
    with(as.list(c(state.seir, parameters)), {
        y.got_infected = y.beta * y.I * y.S / y.N
        y.got_infectious = a * y.E
        y.got_removed = gamma * y.I

        y.deltaS = -y.got_infected
        y.deltaE = y.got_infected - y.got_infectious
        y.deltaI = y.got_infectious - y.got_removed
        y.deltaR = y.got_removed

        o.got_infected = (o.beta * o.I + yo.beta * y.I) * o.S / o.N
        o.got_infectious = a * o.E
        o.got_removed = gamma * o.I

        o.deltaS = -o.got_infected
        o.deltaE = o.got_infected - o.got_infectious
        o.deltaI = o.got_infectious - o.got_removed
        o.deltaR = o.got_removed

        c(y.deltaS, y.deltaE, y.deltaI, y.deltaR,
          o.deltaS, o.deltaE, o.deltaI, o.deltaR)
    })
}

ode.simstep <- function(state, i, y.beta, o.beta, yo.beta, a, gamma)
{
    state.seir <- as.list(state)
    names(state.seir) <- c("y.S", "y.E", "y.I", "y.R",
                           "o.S", "o.E", "o.I", "o.R")

    delta = ode.simstep(i, state.seir, list(y.beta = y.beta,
                                            o.beta = o.beta,
                                            yo.beta = yo.beta,
                                            a = a,
                                            gamma = gamma))

    state + delta
}

simstep <- function(state, y.beta, o.beta, yo.beta, a, gamma)
{
    i = state$i

    deltas <- odesimstepc(state$y.S[i], state$y.E[i], state$y.I[i], state$y.R[i],
                          state$o.S[i], state$o.E[i], state$o.I[i], state$o.R[i],
                          y.beta, o.beta, yo.beta, a, gamma, y.N, o.N)

    state$y.S[i + 1] = state$y.S[i] + deltas[1]
    state$y.E[i + 1] = state$y.E[i] + deltas[2]
    state$y.I[i + 1] = state$y.I[i] + deltas[3]
    state$y.R[i + 1] = state$y.R[i] + deltas[4]

    state$o.S[i + 1] = state$o.S[i] + deltas[5]
    state$o.E[i + 1] = state$o.E[i] + deltas[6]
    state$o.I[i + 1] = state$o.I[i] + deltas[7]
    state$o.R[i + 1] = state$o.R[i] + deltas[8]

    i = i + 1
    state$i = i

    state
}

calcbetas.age <- function(time, betay0, betayt, betao0, betaot, betayo0, betayot)
{
    betai = time - lockdown_offset

    y.beta = o.beta = yo.beta = 0
    
    if (betai < 0) {
        y.beta = betay0
        o.beta = betao0
        yo.beta = betayo0
    } else {
        if (betai >= lockdown_transition_period) {
            y.beta = betayt
            o.beta = betaot
            yo.beta = betayot
        } else {
            fract0 = (lockdown_transition_period - betai) / lockdown_transition_period
            fractt = betai / lockdown_transition_period
            y.beta = fract0 * betay0 + fractt * betayt
            o.beta = fract0 * betao0 + fractt * betaot
            yo.beta = fract0 * betayo0 + fractt * betayot
        }
    }

    c(y.beta, o.beta, yo.beta)
}

calculateModel <- function(params, period)
{
    betay0 <- params[1]
    betayt <- params[2]
    betao0 <- params[3]
    betaot <- params[4]
    betayo0 <- params[5]
    betayot <- params[6]
    hosp_rate <- params[7]
    y.died_rate <- params[8]
    o.died_rate <- params[9]
    y.hosp_latency <- params[10]
    y.died_latency <- params[11]
    o.hosp_latency <- params[12]
    o.died_latency <- params[13]
    Tinf <- params[14]
    Tinc <- params[15]
    mort_lockdown_threshold <- params[16]

    o.hosp_rate_factor = params[17] * o.died_rate / y.died_rate
    y.hosp_rate_factor = 1

    a <- 1 / Tinc
    gamma <- 1 / Tinf

    ## convolution profile to infer hospitalisation count
    y.hosp_cv_profile = calcConvolveProfile(-y.hosp_latency, 5)
    o.hosp_cv_profile = calcConvolveProfile(-o.hosp_latency, 5)

    ## convolution profile to infer dead count
    y.died_cv_profile = calcConvolveProfile(-y.died_latency, 5)
    o.died_cv_profile = calcConvolveProfile(-o.died_latency, 5)

    ## additional padding in state vectors to prevent out-of-bounds
    ## in convolution
    padding = max(-y.hosp_cv_profile$kbegin, -y.died_cv_profile$kbegin,
                  -o.hosp_cv_profile$kbegin, -o.died_cv_profile$kbegin) + 1

    state <- NULL

    state$y.S <- rep(y.N - Initial, padding + period)
    state$y.E <- rep(Initial, padding + period)
    state$y.I <- rep(0, padding + period)
    state$y.R <- rep(0, padding + period)
    state$y.hosp <- rep(0, padding + period)
    state$y.died <- rep(0, padding + period)

    state$o.S <- rep(o.N, padding + period)
    state$o.E <- rep(0, padding + period)
    state$o.I <- rep(0, padding + period)
    state$o.R <- rep(0, padding + period)
    state$o.hosp <- rep(0, padding + period)
    state$o.died <- rep(0, padding + period)

    state$i <- padding + 1

    data_offset = InvalidDataOffset

    y.hr = y.hosp_rate_factor * hosp_rate
    o.hr = o.hosp_rate_factor * hosp_rate

    y.dr = y.died_rate
    o.dr = o.died_rate

    y.beta = o.beta = yo.beta = 0

    if (Tinf > 1 && Tinc > 1) {
        for (i in (padding + 1):(padding + period)) {
            betas = calcbetas.age(i - data_offset,
                                  betay0, betayt,
                                  betao0, betaot,
                                  betayo0, betayot)
            
            state <- simstep(state, betas[1], betas[2], betas[3], a, gamma)

            s = convolute(state$y.S, i, y.hosp_cv_profile)
            state$y.hosp[i] <- (y.N - s) * y.hr

            r = convolute(state$y.R, i, y.died_cv_profile)
            state$y.died[i] <- r * y.dr

            s = convolute(state$o.S, i, o.hosp_cv_profile)
            state$o.hosp[i] <- (o.N - s) * o.hr

            r = convolute(state$o.R, i, o.died_cv_profile)
            state$o.died[i] <- r * o.dr

            if (data_offset == InvalidDataOffset &&
                (state$y.died[i] + state$o.died[i]) > mort_lockdown_threshold) {
                data_offset = i - lockdown_offset
            }
        }
    }

    state$y.deadi <- c(state$y.died[1], diff(state$y.died))
    state$y.hospi <- c(state$y.hosp[1], diff(state$y.hosp))

    state$o.deadi <- c(state$o.died[1], diff(state$o.died))
    state$o.hospi <- c(state$o.hosp[1], diff(state$o.hosp))

    state$padding <- padding
    state$offset <- data_offset

    state
}

transformParams <- function(params)
{
    result = params
    result[8] = exp(params[8] + params[7])
    result[9] = exp(params[9] + params[7])
    result[7] = exp(params[7])
    result[17] = exp(params[17])

    result
}

invTransformParams <- function(posterior)
{
    posterior$HR = exp(posterior$logHR)
    posterior$y.IFR = exp(posterior$logHRyDR + posterior$logHR)
    posterior$o.IFR = exp(posterior$logHRoDR + posterior$logHR)
    posterior$fHo = exp(posterior$logfHo)

    ## Additional quantitites of interest
    posterior$y.R0 = posterior$betay0 * posterior$Tinf
    posterior$y.Rt = posterior$betayt * posterior$Tinf
    posterior$o.R0 = posterior$betao0 * posterior$Tinf
    posterior$o.Rt = posterior$betaot * posterior$Tinf
    posterior$yo.R0 = posterior$betayo0 * posterior$Tinf
    posterior$yo.Rt = posterior$betayot * posterior$Tinf
    posterior$y.Et = posterior$y.Rt / posterior$y.R0
    posterior$o.Et = posterior$o.Rt / posterior$o.R0
    posterior$yo.Et = posterior$yo.Rt / posterior$yo.R0
    posterior$o.H = posterior$HR / posterior$fHo

    posterior
}

calcNominalState <- function(state)
{
    state$hosp <- state$y.hosp + state$o.hosp
    state$died <- state$y.died + state$o.died
    state$hospi <- state$y.hospi + state$o.hospi
    state$deadi <- state$y.deadi + state$o.deadi

    state
}

fit.paramnames <- c("betay0", "betayt",
                    "betao0", "betaot",
                    "betayo0", "betayot",
                    "logHR", "logHRyDR", "logHRoDR",
                    "y.HL", "y.DL", "o.HL", "o.DL",
                    "Tinf", "Tinc",
                    "lockdownmort", "logfHo")

## e.g. used for plotting time series to oversee sampling
keyparamnames <- c("y.R0", "y.Rt", "y.IFR", "o.IFR", "Tinf", "betay0", "betayt")
fitkeyparamnames <- c("betay0", "betayt", "betao0", "logHR", "logHRyDR", "logHRoDR", "Tinf", "Tinc")

## log likelihood function for fitting this model to observed data:
##   y.dhospi, o.dhospi, y.dmorti, o.dmorti
calclogl <- function(params) {
    betay0 <- params[1]
    betayt <- params[2]
    betao0 <- params[3]
    betaot <- params[4]
    betayo0 <- params[5]
    betayot <- params[6]
    hosp_rate <- params[7]
    y.died_rate <- params[8]
    o.died_rate <- params[9]
    y.hosp_latency <- params[10]
    y.died_latency <- params[11]
    o.hosp_latency <- params[12]
    o.died_latency <- params[13]
    Tinf <- params[14]
    Tinc <- params[15]
    mort_lockdown_threshold <- params[16]

    if (betay0 < 1E-10) {
        ## print(paste("invalid betay0", betay0))
        return(-Inf)
    }

    if (betayt < 1E-10) {
        ## print(paste("invalid betayt", betayt))
        return(-Inf)
    }

    if (betao0 < 1E-10) {
        ## print(paste("invalid betao0", betao0))
        return(-Inf)
    }

    if (betaot < 1E-10) {
        ## print(paste("invalid betaot", betaot))
        return(-Inf)
    }

    if (betayo0 < 1E-10) {
        ## print(paste("invalid betayo0", betayo0))
        return(-Inf)
    }

    if (betayot < 1E-10) {
        ## print(paste("invalid betayot", betayot))
        return(-Inf)
    }

    if (y.hosp_latency < 2 || y.hosp_latency > 30) {
        print(paste("invalid y.hosp_latency", y.hosp_latency))
        return(-Inf)
    }

    if (y.died_latency < 2 || y.died_latency > 30) {
        print(paste("invalid y.died_latency", y.died_latency))
        return(-Inf)
    }

    if (o.hosp_latency < -20 || o.hosp_latency > 30) {
        print(paste("invalid o.hosp_latency", o.hosp_latency))
        return(-Inf)
    }

    if (o.died_latency < -20 || o.died_latency > 30) {
        print(paste("invalid o.died_latency", o.died_latency))
        return(-Inf)
    }

    state <<- calculateModel(params, FitTotalPeriod)
    
    if (state$offset == InvalidDataOffset)
        state$offset = 1

    logPriorP <- 0

    ##
    ## Also put here rather uninformative priors to avoid complete flatness
    ##
    logPriorP <- logPriorP + dnorm(betay0, 1, 10, log=T)
    logPriorP <- logPriorP + dnorm(betayt, 1, 10, log=T)
    logPriorP <- logPriorP + dnorm(betao0, 1, 10, log=T)
    logPriorP <- logPriorP + dnorm(betaot, 1, 10, log=T)
    logPriorP <- logPriorP + dnorm(betayo0, 1, 10, log=T)
    logPriorP <- logPriorP + dnorm(betayot, 1, 10, log=T)
    
    ##
    ## Only informative priors are applied
    ##

    ## Priors on IFR:
    ##  IFR can typically not be inferred from the data
    ##  Here we list IFR's based on Verity et al., with different variances

    ## Stronger priors based on Verity
    ## estBetaParams(0.0009, 0.0003^2)
    logPriorP <- logPriorP + dbeta(y.died_rate, 8.991, 9981.009, log=T)
    ## estBetaParams(0.03, 0.01^2)
    logPriorP <- logPriorP + dbeta(o.died_rate, 8.7, 281.3, log=T)

    ## Weaker priors based on Verity
    ## estBetaParams(0.0009, 0.0005^2)
    ##logPriorP <- logPriorP + dbeta(y.died_rate, 3.236, 3592.524, log=T)
    ## estBetaParams(0.03, 0.016^2)
    ##logPriorP <- logPriorP + dbeta(o.died_rate, 3.113, 100.6, log=T)
    
    logPriorP <- logPriorP + dnorm(y.hosp_latency, mean=10, sd=10, log=T)
    logPriorP <- logPriorP + dnorm(y.died_latency, mean=10, sd=10, log=T)

    logPriorP <- logPriorP + dnorm(o.hosp_latency, mean=10, sd=10, log=T)
    logPriorP <- logPriorP + dnorm(o.died_latency, mean=10, sd=10, log=T)

    ##
    ## Based on literature estimates
    ##
    logPriorP <- logPriorP + dnorm(Tinc, mean=5, sd=3, log=T)
    logPriorP <- logPriorP + dnorm(Tinf, mean=5, sd=3, log=T)

    ##
    ## Total deaths at lockdown is also the result of an observation
    ##

    loglLD <- dnbinom(total_deaths_at_lockdown, mu=pmax(0.1, mort_lockdown_threshold),
                      size=mort_nbinom_size, log=T)

    dstart <- state$offset
    dend <- state$offset + length(dhospi) - 1

    if (dend > length(state$y.hospi)) {
        print("=========================== Increase FitTotalPeriod ===================")
        dend <- length(state$y.hospi)
        dstart <- dend - length(dhospi) + 1
    }

    y.loglH <- sum(dnbinom(y.dhospi,
                           mu=pmax(0.1, state$y.hospi[dstart:dend]),
                           size=hosp_nbinom_size, log=T))
    o.loglH <- sum(dnbinom(o.dhospi,
                           mu=pmax(0.1, state$o.hospi[dstart:dend]),
                           size=hosp_nbinom_size, log=T))

    dend <- state$offset + length(y.dmorti) - 1

    if (dend > length(state$y.deadi)) {
       print("=========================== Increase FitTotalPeriod ===================")
       dend <- length(state$y.deadi)
       dstart <- dend - length(y.dmorti) + 1
    }

    y.loglD <- sum(dnbinom(y.dmorti,
                           mu=pmax(0.1, state$y.deadi[dstart:dend]),
                           size=mort_nbinom_size, log=T))

    o.loglD <- sum(dnbinom(o.dmorti,
                           mu=pmax(0.1, state$o.deadi[dstart:dend]),
                           size=mort_nbinom_size, log=T))

    it <<- it + 1

    result <- logPriorP + loglLD + y.loglH + o.loglH + y.loglD + o.loglD
    
    if (it %% 1000 == 0) {
        print(params)
	print(c(it, result))
        state <<- calcNominalState(state)
	graphs()
    }
    
    result
}

init <- c(3.6, 0.7, 3.04, 0.12, 0.27, 0.09,
          log(0.01), log(0.3), log(1),
          14, 13, 14, 13, 1.3, 7, total_deaths_at_lockdown, log(25))
scales <- c(0.15, 0.05, 0.15, 0.05, 0.15, 0.05,
            0.05, 0.05, 0.05,
            1, 1, 1, 1, 0.3, 0.3, total_deaths_at_lockdown / 20, 0.1)
