source("model.R")
require(Rcpp)

cppFunction('NumericVector odesimstepc(double yS, double yE, double yI, double yR, double oS, double oE, double oI, double oR, double ybeta, double obeta, double yobeta, double a, double gamma, int yN, int oN) {
   const int LOOPS = 10;
   const double Ts = 1.0/LOOPS;

   NumericVector out(10);

   for (int l = 0; l < LOOPS; ++l) {
     const double ygot_infected = std::min(1.0, ybeta * yI / yN) * yS;
     const double ygot_infectious = a * yE;
     const double ygot_removed = gamma * yI;

     const double ydeltaS = -ygot_infected;
     const double ydeltaE = ygot_infected - ygot_infectious;
     const double ydeltaI = ygot_infectious - ygot_removed;
     const double ydeltaR = ygot_removed;

     const double ogot_infected = std::min(1.0, obeta * oI / oN + yobeta * yI / yN) * oS;
     const double ogot_infectious = a * oE;
     const double ogot_removed = gamma * oI;

     const double odeltaS = -ogot_infected;
     const double odeltaE = ogot_infected - ogot_infectious;
     const double odeltaI = ogot_infectious - ogot_removed;
     const double odeltaR = ogot_removed;

     if (l == 0) {
        out[8] = 1/gamma * (ygot_infected + ogot_infected) / (yI + oI);
        out[9] = out[8] * (yN + oN) / (yS + oS);
     }

     yS += ydeltaS * Ts;
     yE += ydeltaE * Ts;
     yI += ydeltaI * Ts;
     yR += ydeltaR * Ts;

     oS += odeltaS * Ts;
     oE += odeltaE * Ts;
     oI += odeltaI * Ts;
     oR += odeltaR * Ts;
   }

   out[0] = yS;
   out[1] = yE;
   out[2] = yI;
   out[3] = yR;
   out[4] = oS;
   out[5] = oE;
   out[6] = oI;
   out[7] = oR;

   return out;
}')

simstep <- function(state, y.beta, o.beta, yo.beta, a, gamma)
{
    i = state$i

    newState <- odesimstepc(state$y.S[i], state$y.E[i], state$y.I[i], state$y.R[i],
                            state$o.S[i], state$o.E[i], state$o.I[i], state$o.R[i],
                            y.beta, o.beta, yo.beta, a, gamma, y.N, o.N)

    state$y.S[i + 1] = newState[1]
    state$y.E[i + 1] = newState[2]
    state$y.I[i + 1] = newState[3]
    state$y.R[i + 1] = newState[4]

    state$o.S[i + 1] = newState[5]
    state$o.E[i + 1] = newState[6]
    state$o.I[i + 1] = newState[7]
    state$o.R[i + 1] = newState[8]
    state$Re[i + 1] = newState[9]
    state$Rt[i + 1] = newState[10]

    i = i + 1
    state$i = i

    state
}

calcbetas.age <- function(time, y.beta0, y.betat, o.beta0, o.betat, yo.beta0, yo.betat, Es)
{
    betai = time - lockdown_offset

    y.beta = o.beta = yo.beta = 0
    
    if (betai < 0) {
        y.beta = y.beta0
        o.beta = o.beta0
        yo.beta = yo.beta0
    } else {
        if (betai >= lockdown_transition_period) {
            y.beta = y.betat
            o.beta = o.betat
            yo.beta = yo.betat
            if (length(Es.time) > 0) {
                for (i in length(Es.time):1) {
                    if (time > Es.time[i]) {
                        k <- (i - 1) * 3 + 1
                        y.beta = max(0.001, Es[k] * y.betat + (1 - Es[k]) * y.beta0)
                        o.beta = max(0.001, Es[k + 1] * o.betat + (1 - Es[k + 1]) * o.beta0)
                        yo.beta = max(0.001, Es[k + 2] * yo.betat + (1 - Es[k + 2]) * yo.beta0)
                        break;
                    }
                }
            }
        } else {
            fract0 = (lockdown_transition_period - betai) / lockdown_transition_period
            fractt = betai / lockdown_transition_period
            y.beta = fract0 * y.beta0 + fractt * y.betat
            o.beta = fract0 * o.beta0 + fractt * o.betat
            yo.beta = fract0 * yo.beta0 + fractt * yo.betat
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
    gamma <- params[14]
    a <- params[15]
    mort_lockdown_threshold <- params[16]
    o.hosp_rate_factor = params[17] * o.died_rate / y.died_rate
    Es <- tail(params, n=-17)

    y.hosp_rate_factor = 1

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

    state$Re <- rep(0, padding + period)
    state$Rt <- rep(0, padding + period)
    
    state$i <- padding + 1

    data_offset = InvalidDataOffset

    y.hr = y.hosp_rate_factor * hosp_rate
    o.hr = o.hosp_rate_factor * hosp_rate

    y.dr = y.died_rate
    o.dr = o.died_rate

    y.beta = o.beta = yo.beta = 0

    for (i in (padding + 1):(padding + period)) {
        betas = calcbetas.age(i - data_offset,
                              betay0, betayt,
                              betao0, betaot,
                              betayo0, betayot,
                              Es)

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
    result[1:6] = exp(params[1:6])
    result[8] = exp(params[8] + params[7])
    result[9] = exp(params[9] + params[7])
    result[14:15] = exp(params[14:15])
    result[7] = exp(params[7])
    result[17] = exp(params[17])

    result
}

invTransformParams <- function(posterior)
{
    posterior$betay0 = exp(posterior$logBetay0)
    posterior$betayt = exp(posterior$logBetayt)
    posterior$betao0 = exp(posterior$logBetao0)
    posterior$betaot = exp(posterior$logBetaot)
    posterior$betayo0 = exp(posterior$logBetayo0)
    posterior$betayot = exp(posterior$logBetayot)
    posterior$HR = exp(posterior$logHR)
    posterior$y.IFR = exp(posterior$logHRyDR + posterior$logHR)
    posterior$o.IFR = exp(posterior$logHRoDR + posterior$logHR)
    posterior$fHo = exp(posterior$logfHo)
    posterior$a = exp(posterior$loga)
    posterior$Tinc = 1/posterior$a
    posterior$gamma = exp(posterior$logGamma)
    posterior$Tinf = 1/posterior$gamma

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
    gamma <- params[14]
    a <- params[15]
    mort_lockdown_threshold <- params[16]
    o.hosp_rate_factor = params[17]
    Es <- tail(params, n=-17)

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
    ## Based on literature estimates, tightened now to converge more quickly
    ##
    Tinc = 1/a
    logPriorP <- logPriorP + dnorm(Tinc, mean=2.5, sd=0.5, log=T)
    Tinf = 1/gamma
    logPriorP <- logPriorP + dnorm(Tinf, mean=2.5, sd=0.5, log=T)

    ## Es: prior is 1 for first one; prior for Es[i] = estimate for Es[i-1]
    ## (thus assume for each subsequent period no change relative to the previous period)
    if (length(Es) > 0) {
        logPriorP <- logPriorP + dnorm(Es[1], mean=1.0, sd=0.1, log=T)
        if (length(Es) > 1) {
            for (eI in 2:length(Es)) {
                logPriorP <- logPriorP + dnorm(Es[eI], mean=Es[eI-1], sd=0.1, log=T)
            }
        }
    }

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

fit.paramnames <- c("logBetay0", "logBetayt",
                    "logBetao0", "logBetaot",
                    "logBetayo0", "logBetayot",
                    "logHR", "logHRyDR", "logHRoDR",
                    "y.HL", "y.DL", "o.HL", "o.DL",
                    "logGamma", "loga",
                    "lockdownmort", "logfHo")

## e.g. used for plotting time series to oversee sampling
keyparamnames <- c("y.R0", "y.Rt", "y.IFR", "o.IFR", "Tinf", "betay0", "betayt")
fitkeyparamnames <- c("logBetay0", "logBetayt", "logBetao0", "logHR", "logHRyDR", "logHRoDR", "logGamma", "loga")

init <- c(log(3.6), log(0.7), log(3.04), log(0.12), log(0.27), log(0.09),
          log(0.01), log(0.3), log(1),
          14, 13, 14, 13, 1, 0.3, total_deaths_at_lockdown, 0.2,
          rep(0.9, length(Es.time) * 3))
scales <- c(1, 1, 1, 1, 1, 1,
            0.05, 0.05, 0.05,
            1, 1, 1, 1, 1, 1, total_deaths_at_lockdown / 20, 0.1,
            rep(0.05, length(Es.time) * 3))

if (length(Es.time) > 0) {
    for (i in 1:length(Es.time)) {
        fit.paramnames <- c(fit.paramnames, paste("y.E", i, sep=""),
                            paste("o.E", i, sep=""),
                            paste("yo.E", i, sep=""))
    }
}
