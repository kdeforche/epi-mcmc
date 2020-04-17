source("model.R")

simstep <- function(state, y.beta, o.beta, yo.beta, a, gamma)
{
    i = state$i

    y.got_infected = y.beta * state$y.I[i] * state$y.S[i] / y.N
    y.got_infectious = a * state$y.E[i]
    y.got_removed = gamma * state$y.I[i]

    y.deltaS = -y.got_infected
    y.deltaE = y.got_infected - y.got_infectious
    y.deltaI = y.got_infectious - y.got_removed
    y.deltaR = y.got_removed

    state$y.S[i + 1] = state$y.S[i] + y.deltaS
    state$y.E[i + 1] = state$y.E[i] + y.deltaE
    state$y.I[i + 1] = state$y.I[i] + y.deltaI
    state$y.R[i + 1] = state$y.R[i] + y.deltaR

    o.got_infected = (o.beta * state$o.I[i] + yo.beta * state$y.I[i]) * state$o.S[i] / o.N
    o.got_infectious = a * state$o.E[i]
    o.got_removed = gamma * state$o.I[i]

    o.deltaS = -o.got_infected
    o.deltaE = o.got_infected - o.got_infectious
    o.deltaI = o.got_infectious - o.got_removed
    o.deltaR = o.got_removed

    state$o.S[i + 1] = state$o.S[i] + o.deltaS
    state$o.E[i + 1] = state$o.E[i] + o.deltaE
    state$o.I[i + 1] = state$o.I[i] + o.deltaI
    state$o.R[i + 1] = state$o.R[i] + o.deltaR

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
    died_rate <- params[8]
    y.hosp_latency <- params[9]
    y.died_latency <- params[10]
    o.hosp_latency <- params[11]
    o.died_latency <- params[12]
    Tinf <- params[13]
    Tinc <- params[14]
    mort_lockdown_threshold <- params[15]
    o.hosp_rate_factor = params[16]

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
    y.dr = y.died_rate_factor * died_rate

    o.hr = o.hosp_rate_factor * hosp_rate
    o.dr = o.died_rate_factor * died_rate

    y.beta = o.beta = yo.beta = 0
    
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

    state$y.deadi <- numeric(length(state$y.died))
    state$y.hospi <- numeric(length(state$y.hosp))

    for (i in (padding + 1):(padding + period)) {
    	state$y.deadi[i] <- if (i == 1) state$y.died[1]
                            else state$y.died[i] - state$y.died[i - 1]
	state$y.hospi[i] <- if (i == 1) state$y.hosp[1]
                            else state$y.hosp[i] - state$y.hosp[i - 1]
    }

    state$o.deadi <- numeric(length(state$o.died))
    state$o.hospi <- numeric(length(state$o.hosp))

    for (i in (padding + 1):(padding + period)) {
    	state$o.deadi[i] <- if (i == 1) state$o.died[1]
                            else state$o.died[i] - state$o.died[i - 1]
	state$o.hospi[i] <- if (i == 1) state$o.hosp[1]
                            else state$o.hosp[i] - state$o.hosp[i - 1]
    }

    state$padding <- padding
    state$offset <- data_offset

    state
}

transformParams <- function(params)
{
    result = params
    result[8] = exp(params[7] + params[8])
    result[7] = exp(params[7])
    result[16] = exp(params[16])

    result
}

invTransformParams <- function(posterior)
{
    posterior$IFR = exp(posterior$logHR + posterior$logHRDR)
    posterior$HR = exp(posterior$logHR)
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
                    "logHR", "logHRDR",
                    "y.HL", "y.DL", "o.HL", "o.DL",
                    "Tinf", "Tinc",
                    "lockdownmort", "logfHo")

## e.g. used for plotting time series to oversee sampling
keyparamnames <- c("y.R0", "y.Rt", "IFR", "Tinf", "betay0", "betayt")

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
    died_rate <- params[8]
    y.hosp_latency <- params[9]
    y.died_latency <- params[10]
    o.hosp_latency <- params[11]
    o.died_latency <- params[12]
    Tinf <- params[13]
    Tinc <- params[14]
    mort_lockdown_threshold <- params[15]

    if (betay0 < 1E-10) {
        return(-Inf)
    }

    if (betayt < 1E-10) {
        return(-Inf)
    }

    if (betao0 < 1E-10) {
        return(-Inf)
    }

    if (betaot < 1E-10) {
        return(-Inf)
    }

    if (betayo0 < 1E-10) {
        return(-Inf)
    }

    if (betayot < 1E-10) {
        return(-Inf)
    }

    if (hosp_rate < 0.0001) {
        print(paste("invalid hosp_rate", hosp_rate))
        return(-Inf)
    }

    if (died_rate < 0.0001) {
        print(paste("invalid died_rate", died_rate))
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

    if (o.hosp_latency < -2 || o.hosp_latency > 30) {
        print(paste("invalid o.hosp_latency", o.hosp_latency))
        return(-Inf)
    }

    if (o.died_latency < -2 || o.died_latency > 30) {
        print(paste("invalid o.died_latency", o.died_latency))
        return(-Inf)
    }

    if (Tinf < 1.1 || Tinf > 20) {
        #print(paste("invalid Tinf", Tinf))
        return(-Inf)
    }

    if (Tinc < 1.1 || Tinc > 20) {
        print(paste("invalid Tinc", Tinc))
        return(-Inf)
    }

    state <<- calculateModel(params, FitTotalPeriod)
    
    if (state$offset == InvalidDataOffset)
       return(-Inf)

    logPriorP <- 0

    ##
    ## Only informative priors are applied
    ##

    ## Priors on IFR:
    ##  IFR can typically not be inferred from the data
    ##  Here we list IFR's based on Verity et al., with different variances

    ## estBetaParams(0.0066, 0.002^2)
    ## logPriorP <- logPriorP + dbeta(died_rate, 10.8, 1627, log=T)

    ## estBetaParams(0.0066, 0.003^2)
    ## logPriorP <- logPriorP + dbeta(died_rate, 4.8, 722.69, log=T)

    ## estBetaParams(0.0066, 0.004^2)
    ## logPriorP <- logPriorP + dbeta(died_rate, 2.697931, 406.0796, log=T)

    ## estBetaParams(0.0066, 0.005^2)
    logPriorP <- logPriorP + dbeta(died_rate, 1.7243, 259.53, log=T)

    ##
    ## These two may not be needed
    ##
    logPriorP <- logPriorP + dnorm(y.hosp_latency, mean=10, sd=20, log=T)
    logPriorP <- logPriorP + dnorm(y.died_latency, mean=10, sd=20, log=T)

    logPriorP <- logPriorP + dnorm(o.hosp_latency, mean=10, sd=20, log=T)
    logPriorP <- logPriorP + dnorm(o.died_latency, mean=10, sd=20, log=T)

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
       return(-Inf)
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
       return(-Inf)
    }

    y.loglD <- sum(dnbinom(y.dmorti,
                           mu=pmax(0.1, state$y.deadi[dstart:dend] /
                                        death_underreporting_factor),
                           size=mort_nbinom_size, log=T))

    o.loglD <- sum(dnbinom(o.dmorti,
                           mu=pmax(0.1, state$o.deadi[dstart:dend] /
                                        death_underreporting_factor),
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
