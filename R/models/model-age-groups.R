source("model.R")

y.N <- 8.55E6
o.N <- N - y.N

## Group-specific factors for dead_rate with respect to estimated died_rate
y.died_rate_factor = 0.09412495/0.66
o.died_rate_factor = 2.34371716/0.66

## Group-specific factors for hosp_rate with respect to estimated hosp_rate
y.hosp_rate_factor = 0.09412495/0.66
o.hosp_rate_factor = 2.34371716/0.66

dmort <- y.dmort + o.dmort
dmorti <- y.dmorti + o.dmorti

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
    hosp_latency <- params[9]
    died_latency <- params[10]
    Tinf <- params[11]
    Tinc <- params[12]
    mort_lockdown_threshold <- params[13]
    ##hosp_rate_change <- params[14]

    a <- 1 / Tinc
    gamma <- 1 / Tinf

    hosp_cv_profile = calcConvolveProfile(-hosp_latency, 5)

    died_cv_profile = calcConvolveProfile(-died_latency, 5)

    padding = floor(max(hosp_latency, died_latency) * 3)

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

    beta_transition = 7
    data_offset = InvalidDataOffset

    y.hr = y.hosp_rate_factor * hosp_rate
    y.dr = y.died_rate_factor * died_rate

    o.hr = o.hosp_rate_factor * hosp_rate
    o.dr = o.died_rate_factor * died_rate

    y.beta = o.beta = yo.beta = 0
    
    for (i in (padding + 1):(padding + period)) {
        betai = i - data_offset - lockdown_offset

        if (betai < 0) {
            y.beta = betay0
            o.beta = betao0
            yo.beta = betayo0
	} else {
            if (betai >= beta_transition) {
                y.beta = betayt
                o.beta = betaot
                yo.beta = betayot
            } else {
                fract0 = (beta_transition - betai) / beta_transition
                fractt = betai / beta_transition
                y.beta = fract0 * betay0 + fractt * betayt
                o.beta = fract0 * betao0 + fractt * betaot
                yo.beta = fract0 * betayo0 + fractt * betayot
            }
	}

	state <- simstep(state, y.beta, o.beta, yo.beta, a, gamma)

    	s = convolute(state$y.S, i, hosp_cv_profile)
	state$y.hosp[i] <- (y.N - s) * y.hr

	r = convolute(state$y.R, i, died_cv_profile)
	state$y.died[i] <- r * y.dr

    	s = convolute(state$o.S, i, hosp_cv_profile)
	state$o.hosp[i] <- (o.N - s) * o.hr

	r = convolute(state$o.R, i, died_cv_profile)
	state$o.died[i] <- r * o.dr

        ## assuming lockdown decision was based on a cumulative mort count, with some
        ## uncertainty on the exact value due to observed cumulative mort count being a
        ## sample
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

    ## cumDiff <- 0
    ## inv_hosp_rate_change <- 1 / hosp_rate_change
    ## for (i in (round(data_offset + 12)):length(state$y.hospi)) {
    ## 	nhospi <- state$y.hospi[i] * inv_hosp_rate_change
    ##     diff <- nhospi - state$y.hospi[i]
    ## 	state$y.hospi[i] <- nhospi
    ##     cumDiff <- cumDiff + diff
    ##     state$y.hosp[i] <- state$y.hosp[i] + cumDiff
    ## }

    ## for (i in (round(data_offset + 12)):length(state$o.hospi)) {
    ## 	nhospi <- state$o.hospi[i] * inv_hosp_rate_change
    ##     diff <- nhospi - state$o.hospi[i]
    ## 	state$o.hospi[i] <- nhospi
    ##     cumDiff <- cumDiff + diff
    ##     state$o.hosp[i] <- state$o.hosp[i] + cumDiff
    ## }

    state$padding <- padding
    state$offset <- data_offset

    state
}

transformParams <- function(params)
{
    result = params
    result[8] = exp(params[7] + params[8])
    result[7] = exp(params[7])

    result
}

fit.paramnames <- c("betay0", "betayt",
                    "betao0", "betaot",
                    "betayo0", "betayot",
                    "logHR", "logHRDR", "HL", "DL", "Tinf", "Tinc",
                    "lockdownmort")

## log likelihood function for fitting this model to observed data:
##   dhospi,
##   y.dmorti, o.dmorti
calclogl <- function(params) {
    betay0 <- params[1]
    betayt <- params[2]
    betao0 <- params[3]
    betaot <- params[4]
    betayo0 <- params[5]
    betayot <- params[6]
    hosp_rate <- params[7]
    died_rate <- params[8]
    hosp_latency <- params[9]
    died_latency <- params[10]
    Tinf <- params[11]

    Tinc <- params[12]
    mort_lockdown_threshold <- params[13]
    ##hosp_rate_change <- params[14]

    if (betay0 < 1E-10) {
        #print(paste("invalid betay0", betay0))
        return(-Inf)
    }

    if (betayt < 1E-10) {
        #print(paste("invalid betay0", betayt))
        return(-Inf)
    }

    if (betao0 < 1E-10) {
        #print(paste("invalid betao0", betao0))
        return(-Inf)
    }

    if (betaot < 1E-10) {
        #print(paste("invalid betaot", betaot))
        return(-Inf)
    }

    if (betayo0 < 1E-10) {
        #print(paste("invalid betayo0", betayo0))
        return(-Inf)
    }

    if (betayot < 1E-10) {
        #print(paste("invalid betayot", betayot))
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

    if (hosp_latency < 2 || hosp_latency > 30) {
        print(paste("invalid hosp_latency", hosp_latency))
        return(-Inf)
    }

    if (died_latency < 2 || died_latency > 30) {
        ## print(paste("invalid died_latency", died_latency))
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

    ##if (hosp_rate_change < 1.0) {
    ##    print(paste("invalid hosp_rate_change", hosp_rate_change))
    ##    return(-Inf)
    ##}

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
    
    logl <- logl + dbeta(died_rate, 10.8, 1627, log=T) # estBetaParams(0.0066, 0.002^2)
    #logl <- logl + dbeta(died_rate, 4.8, 722.69, log=T) # estBetaParams(0.0066, 0.003^2)
    #logl <- logl + dbeta(died_rate, 2.697931, 406.0796, log=T) # estBetaParams(0.0066, 0.004^2)
    #logl <- logl + dbeta(died_rate, 1.7243, 259.53, log=T) # estBetaParams(0.0066, 0.005^2)
    #logl <- logl + dbeta(died_rate, 0.9, 0.9, log=T) # https://stats.stackexchange.com/questions/297901/choosing-between-uninformative-beta-priors
    logl <- logl + dnorm(hosp_latency, mean=10, sd=20, log=T)
    logl <- logl + dnorm(died_latency, mean=10, sd=20, log=T)
    logl <- logl + dnorm(Tinc, mean=5, sd=3, log=T)
    logl <- logl + dnorm(Tinf, mean=5, sd=3, log=T)
    logl <- logl + dnbinom(total_deaths_at_lockdown, mu=pmax(0.1, mort_lockdown_threshold), size=mort_nbinom_size, log=T)
    #logl <- logl + dnorm(hosp_rate_change, mean=1, sd=0.3, log=T)

    dstart <- state$offset
    dend <- state$offset + length(dhospi) - 1

    if (dend > length(state$y.hospi)) {
       print("=========================== Increase FitTotalPeriod ===================")
       return(-Inf)
    }

    loglA <- sum(dnbinom(dhospi,
                         mu=pmax(0.1, state$y.hospi[dstart:dend] +
                                      state$o.hospi[dstart:dend]),
                         size=hosp_nbinom_size, log=T))

    dend <- state$offset + length(y.dmorti) - 1

    if (dend > length(state$y.deadi)) {
       print("=========================== Increase FitTotalPeriod ===================")
       return(-Inf)
    }

    y.loglB <- sum(dnbinom(y.dmorti,
                           mu=pmax(0.1, state$y.deadi[dstart:dend] / death_underreporting_factor),
                           size=mort_nbinom_size, log=T))

    o.loglB <- sum(dnbinom(o.dmorti,
                           mu=pmax(0.1, state$o.deadi[dstart:dend] / death_underreporting_factor),
                           size=mort_nbinom_size, log=T))

    it <<- it + 1
    if (it %% 1000 == 0) {
        print(params)
	print(c(it, logl + loglA + y.loglB + o.loglB))
        state$hosp <<- state$y.hosp + state$o.hosp
        state$died <<- state$y.died + state$o.died
        state$hospi <<- state$y.hospi + state$o.hospi
        state$deadi <<- state$y.deadi + state$o.deadi
	graphs()
    }

    logl + loglA + y.loglB + o.loglB
}
