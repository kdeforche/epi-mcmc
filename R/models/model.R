InvalidDataOffset <- 10000
Initial <- 1

calcConvolveProfile <- function(latency, latency_sd)
{
    kend = min(0, ceiling(latency + latency_sd * 1.5))
    kbegin = min(kend - 1, floor(latency - latency_sd * 1.5))

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

simstep <- function(state, N, Ne, beta, a, gamma)
{
    i = state$i

    EIR = N - state$S[i]

    got_infected = beta * min(1, state$I[i] / Ne) * max(0, Ne - EIR)
    got_infectious = a * state$E[i]
    got_removed = gamma * state$I[i]

    deltaS = -got_infected
    deltaE = got_infected - got_infectious
    deltaI = got_infectious - got_removed
    deltaR = got_removed

    state$Re[i + 1] = got_infected / state$I[i] / gamma
    state$Rt[i + 1] = beta / gamma
    state$S[i + 1] = state$S[i] + deltaS
    state$E[i + 1] = state$E[i] + deltaE
    state$I[i + 1] = state$I[i] + deltaI
    state$R[i + 1] = state$R[i] + deltaR

    i = i + 1
    state$i = i

    state
}

calcbeta <- function(time, beta0, betat)
{
    betai = time - lockdown_offset

    beta = 0
    
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

    beta
}

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
    Nef <- params[10]

    Ne <- N * Nef
    a <- 1 / Tinc
    gamma <- 1 / Tinf

    hosp_cv_profile = calcConvolveProfile(-hosp_latency, 5)
    died_cv_profile = calcConvolveProfile(-died_latency, 5)

    padding = max(-hosp_cv_profile$kbegin, -died_cv_profile$kbegin) + 1

    state <- NULL
    state$Re <- rep(0, padding + period)
    state$Rt <- rep(0, padding + period)
    state$S <- rep(N-Initial, padding + period)
    state$E <- rep(Initial, padding + period)
    state$I <- rep(0, padding + period)
    state$R <- rep(0, padding + period)
    state$hosp <- rep(0, padding + period)
    state$died <- rep(0, padding + period)
    state$i <- padding + 1

    data_offset = InvalidDataOffset

    if (Tinf > 1 && Tinc > 1) {
        beta = beta0	
        for (i in (padding + 1):(padding + period)) {
            beta = calcbeta(i - data_offset, beta0, betat)

            state <- simstep(state, N, Ne, beta, a, gamma)

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
    }

    state$deadi <- c(state$died[1], diff(state$died))
    state$hospi <- c(state$hosp[1], diff(state$hosp))

    state$padding <- padding
    state$offset <- data_offset

    state
}

calcNominalState <- function(state)
{
    state
}

transformParams <- function(params)
{
    result = params
    result[4] = exp(params[3] + params[4])
    result[3] = exp(params[3])
    result[10] = 1 ## Nef

    result
}

invTransformParams <- function(posterior)
{
    posterior$IFR = exp(posterior$logHR + posterior$logHRDR)
    posterior$HR = exp(posterior$logHR)

    ## Additional quantitites of interest
    posterior$R0 = posterior$beta0 * posterior$Tinf
    posterior$Rt = posterior$betat * posterior$Tinf
    posterior$Et = posterior$Rt / posterior$R0

    posterior
}

fit.paramnames <- c("beta0", "betat", "logHR", "logHRDR", "HL", "DL",
                    "Tinf", "Tinc", "lockdownmort")
keyparamnames <- c("R0", "Rt", "IFR", "Tinf", "Tinc", "beta0", "betat")
fitkeyparamnames <- c("beta0", "betat", "logHR", "logHRDR", "Tinf", "Tinc")

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
    Nef <- params[10]

    if (Nef < 0 || Nef > 1) {
        return(-Inf)
    }

    if (beta0 < 1E-10) {
        ## print(paste("invalid beta0", beta0))
        return(-Inf)
    }

    if (betat < 1E-10) {
        ## print(paste("invalid betat", betat))
        return(-Inf)
    }

    if (hosp_latency < 0 || hosp_latency > 30) {
        ##print(paste("invalid hosp_latency", hosp_latency))
        return(-Inf)
    }

    if (died_latency < 0 || died_latency > 30) {
        ##print(paste("invalid died_latency", died_latency))
        return(-Inf)
    }

    state <<- calculateModel(params, FitTotalPeriod)

    if (state$offset == InvalidDataOffset)
        state$offset = 1

    logPriorP <- 0

    logPriorP <- logPriorP + dnorm(beta0, 1, 10, log=T)
    logPriorP <- logPriorP + dnorm(betat, 1, 10, log=T)

    ## estBetaParams(0.0066, 0.005^2)
    ## logPriorP <- logPriorP + dbeta(died_rate, 1.7243, 259.53, log=T)
    ## estBetaParams(0.0066, 0.002^2)
    logPriorP <- logPriorP + dbeta(died_rate, 10.811, 1627.298, log=T)

    logPriorP <- logPriorP + dnorm(hosp_latency, mean=10, sd=10, log=T)
    logPriorP <- logPriorP + dnorm(died_latency, mean=10, sd=10, log=T)

    logPriorP <- logPriorP + dnorm(Tinc, mean=5, sd=3, log=T)
    logPriorP <- logPriorP + dnorm(Tinf, mean=5, sd=3, log=T)

    loglLD <- dnbinom(total_deaths_at_lockdown, mu=pmax(0.1, mort_lockdown_threshold),
                      size=mort_nbinom_size, log=T)

    dstart <- state$offset
    dend <- state$offset + length(dhospi) - 1

    if (dend > length(state$hospi)) {
        print("=========================== Increase FitTotalPeriod ===================")
        dend <- length(state$hospi)
        dstart <- dend - length(dhospi) + 1
    }

    loglH <- sum(dnbinom(dhospi,
                         mu=pmax(0.1, state$hospi[dstart:dend]),
                         size=hosp_nbinom_size, log=T))

    dend <- state$offset + length(dmorti) - 1

    if (dend > length(state$deadi)) {
        print("=========================== Increase FitTotalPeriod ===================")
        dend <- length(state$deadi)
        dstart <- dend - length(dmorti) + 1
    }

    loglD <- sum(dnbinom(dmorti,
                         mu=pmax(0.1, state$deadi[dstart:dend]), size=mort_nbinom_size, log=T))

    it <<- it + 1

    result <- logPriorP + loglLD + loglH + loglD
    
    if (it %% 1000 == 0) {
        print(params)
	print(c(it, result))
	graphs()
    }

    result
}

init <- c(2.3, 0.5, log(0.02), log(0.3), 10, 9, 2, 5, total_deaths_at_lockdown)
scales <- c(0.15, 0.05, 0.05, 0.1, 1, 1, 0.3, 0.3, total_deaths_at_lockdown / 20)
