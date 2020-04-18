InvalidDataOffset <- 10000
Initial <- 1

calcConvolveProfile <- function(latency, latency_sd)
{
    kbegin = floor(latency - latency_sd * 1.5)
    kend = min(0, ceiling(latency + latency_sd * 1.5))

    result = NULL
    result$kbegin = kbegin
    result$kend = kend
    result$values = numeric(kend - kbegin)
    i = 1
    for (k in kbegin:kend) {
        result$values[i] = pnorm(k-0.5, mean=latency, sd=latency_sd) -
            pnorm(k+0.5, mean=latency, sd=latency_sd)
        i = i + 1
    }

    result$values = result$values / sum(result$values)    

    result
}    

convolute <- function(values, i, profile)
{
    i1 <- i + profile$kbegin
    i2 <- i + profile$kend

    (profile$values %*% values[i1:i2])[1,1]
}

simstep <- function(state, deltaT, N, beta, a, gamma)
{
    i = state$i
    got_infected = beta * state$I[i] / N * state$S[i]
    got_infectious = a * state$E[i]
    got_removed = gamma * state$I[i]

    deltaS = -got_infected
    deltaE = got_infected - got_infectious
    deltaI = got_infectious - got_removed
    deltaR = got_removed

    state$S[i + 1] = state$S[i] + deltaS * deltaT
    state$E[i + 1] = state$E[i] + deltaE * deltaT
    state$I[i + 1] = state$I[i] + deltaI * deltaT
    state$R[i + 1] = state$R[i] + deltaR * deltaT

    i = i + 1
    state$i = i

    state
}

model.paramnames <- c("beta0", "betat", "HR", "DR", "HL", "DL", "Tinf", "Tinc", "lockdownmort", "WZC")

calculateModel <- function(params, period)
{
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

    a <- 1 / Tinc
    gamma <- 1 / Tinf

    hosp_cv_profile = calcConvolveProfile(-hosp_latency, 5)
    died_cv_profile = calcConvolveProfile(-died_latency, 5)

    padding = max(-hosp_cv_profile$kbegin, -died_cv_profile$kbegin) + 1

    state <- NULL
    state$S <- rep(N-Initial, padding + period)
    state$E <- rep(Initial, padding + period)
    state$I <- rep(0, padding + period)
    state$R <- rep(0, padding + period)
    state$hosp <- rep(0, padding + period)
    state$died <- rep(0, padding + period)
    state$i <- padding + 1

    data_offset = InvalidDataOffset

    beta = beta0	
    for (i in (padding + 1):(padding + period)) {
        betai = i - data_offset - lockdown_offset
	if (betai < 0) {
            beta = beta0
	} else {
            if (betai >= lockdown_transition_period) {
                beta = betat
            } else {
                fract0 = (lockdown_transition_period - betai) / lockdown_transition_period
                fractt = betai / lockdown_transition_period
                beta = fract0 * beta0 + fractt * betat
            }
	}

	state <- simstep(state, 1, N, beta, a, gamma)

    	s = convolute(state$S, i, hosp_cv_profile)
	state$hosp[i] <- (N - s) * hosp_rate

	r = convolute(state$R, i, died_cv_profile)
	state$died[i] <- r * died_rate

        ## assuming lockdown decision was based on a cumulative mort count, with some
        ## uncertainty on the exact value due to observed cumulative mort count being a
        ## sample
	if (data_offset == InvalidDataOffset && state$died[i] > mort_lockdown_threshold) {
            data_offset = i - lockdown_offset
	}
    }

    state$deadi <- numeric(length(state$died))
    state$hospi <- numeric(length(state$hosp))

    for (i in (padding + 1):(padding + period)) {
    	state$deadi[i] <- if (i == 1) state$died[1] else state$died[i] - state$died[i - 1]
	state$hospi[i] <- if (i == 1) state$hosp[1] else state$hosp[i] - state$hosp[i - 1]
    }

    ## FIXME only applies to Belgian data set -- hospitals reduced the rate at which
    ## they were accepting elderly people on March 22
    cumDiff <- 0
    inv_hosp_rate_change <- 1 / hosp_rate_change
    for (i in (round(data_offset + 12)):length(state$hospi)) {
    	nhospi <- state$hospi[i] * inv_hosp_rate_change
	diff <- nhospi - state$hospi[i]
    	state$hospi[i] <- nhospi
	cumDiff <- cumDiff + diff
	state$hosp[i] <- state$hosp[i] + cumDiff
    }

    state$padding <- padding
    state$offset <- data_offset

    state
}

calcNominalState <- function(state)
{
    state
}

fit.paramnames <- c("beta0", "betat", "logHR", "logHRDR", "HL", "DL", "Tinf", "Tinc",
                    "lockdownmort", "logWZC")

## e.g. used for plotting time series to oversee sampling
keyparamnames <- c("R0", "Rt", "IFR", "Tinf", "beta0", "betat")

transformParams <- function(params)
{
    result = params
    result[4] = exp(params[3] + params[4])
    result[3] = exp(params[3])
    result[10] = exp(params[10])

    result
}

invTransformParams <- function(posterior)
{
    posterior$IFR = exp(posterior$logHR + posterior$logHRDR)
    posterior$HR = exp(posterior$logHR)
    posterior$WZC = exp(posterior$logWZC)

    ## Additional quantitites of interest
    posterior$R0 = posterior$beta0 * posterior$Tinf
    posterior$Rt = posterior$betat * posterior$Tinf
    posterior$Et = posterior$Rt / posterior$R0

    posterior
}


## log likelihood function for fitting this model to observed data:
##   dhospi and dmorti
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
