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

simstep <- function(state, N, beta, a, gamma)
{
    i = state$i

    got_infected = beta * state$I[i] / N * state$S[i]
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

##
## Beta is a function of time:
##  t < lockdown_offset : beta0
##  t > lockdown_offset + lockdown_transition_period : betat
##  t > Es.time[i] : Es[i] * betat + (1 - Es[i]) * beta0

if (!exists("Es.time")) {
    Es.time <- c()
}

calcpar <- function(time, par0, part, Es)
{
    pari = time - lockdown_offset

    par = 0
    
    if (pari < 0) {
        par = par0
    } else {
        if (pari >= lockdown_transition_period) {
            par = part
            if (length(Es.time) > 0) {
                for (i in length(Es.time):1) {
                    if (time > Es.time[i]) {
                        par = max(0.001, Es[i] * part + (1 - Es[i]) * par0)
                        break;
                    }
                }
            }
        } else {
            fract0 = (lockdown_transition_period - pari) / lockdown_transition_period
            fractt = pari / lockdown_transition_period
            par = fract0 * par0 + fractt * part
        }
    }

    par
}

calculateModel <- function(params, period)
{
    beta <- params[1]
    hosp_rate <- params[2]
    died_rate <- params[3]
    hosp_latency <- params[4]
    died_latency <- params[5]
    Tinf0 <- params[6]
    Tinft <- params[7]
    Tinc <- params[8]
    mort_lockdown_threshold <- params[9]
    Es <- tail(params, n=-9)

    a <- 1 / Tinc

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

    if (Tinf0 > 1 && Tinft > 1 && Tinc > 1) {
        for (i in (padding + 1):(padding + period)) {
            Tinf = calcpar(i - data_offset, Tinf0, Tinft, Es)
            gamma <- 1 / Tinf

            state <- simstep(state, N, beta, a, gamma)

            s = convolute(state$S, i, hosp_cv_profile)
            state$hosp[i] <- (N - s) * hosp_rate

            ## died:
            ##  from those that became infectious
            
            r = convolute(state$I + state$R, i, died_cv_profile)
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
    result[3] = exp(params[2] + params[3])
    result[2] = exp(params[2])

    result
}

invTransformParams <- function(posterior)
{
    posterior$IFR = exp(posterior$logHR + posterior$logHRDR)
    posterior$HR = exp(posterior$logHR)

    ## Additional quantitites of interest
    posterior$R0 = posterior$beta * posterior$Tinf0
    posterior$Rt = posterior$beta * posterior$Tinft

    posterior
}

## log likelihood function for fitting this model to observed data:
##   dhospi and dmorti
calclogl <- function(params) {
    beta <- params[1]
    hosp_rate <- params[2]
    died_rate <- params[3]
    hosp_latency <- params[4]
    died_latency <- params[5]
    Tinf0 <- params[6]
    Tinft <- params[7]
    Tinc <- params[8]
    mort_lockdown_threshold <- params[9]
    Es <- tail(params, n=-9)

    if (beta < 1E-10) {
        ## print(paste("invalid beta0", beta0))
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

    logPriorP <- logPriorP + dexp(beta, 0.01, log=T)

    ## estBetaParams(0.0066, 0.005^2)
    ## logPriorP <- logPriorP + dbeta(died_rate, 1.7243, 259.53, log=T)
    ## estBetaParams(0.0066, 0.002^2)
    ## logPriorP <- logPriorP + dbeta(died_rate, 10.811, 1627.298, log=T)
    ## estBetaParams(0.0066, 0.001^2)
    logPriorP <- logPriorP + dbeta(died_rate, 43.2659, 6512.174, log=T)

    logPriorP <- logPriorP + dnorm(hosp_latency, mean=10, sd=10, log=T)
    logPriorP <- logPriorP + dnorm(died_latency, mean=10, sd=10, log=T)

    logPriorP <- logPriorP + dnorm(Tinc, mean=5, sd=3, log=T)
    logPriorP <- logPriorP + dnorm(Tinf0, mean=5, sd=5, log=T)
    logPriorP <- logPriorP + dnorm(Tinft, mean=5, sd=5, log=T)

    for (e in Es) {
        logPriorP <- logPriorP + dnorm(e, mean=0.9, sd=0.1, log=T)
    }

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

fit.paramnames <- c("beta", "logHR", "logHRDR", "HL", "DL",
                    "Tinf0", "Tinft", "Tinc", "lockdownmort")
keyparamnames <- c("beta", "IFR", "Tinf0", "Tinf1", "Tinc", "R0", "Rt")
fitkeyparamnames <- c("beta", "logHR", "logHRDR", "Tinf0", "Tinft", "Tinc")

init <- c(1, log(0.02), log(0.3), 10, 9, 5, 5, 5, total_deaths_at_lockdown,
          rep(0.9, length(Es.time)))
scales <- c(1, 0.05, 0.1, 1, 1, 0.3, 0.3, 0.3, total_deaths_at_lockdown / 20,
            rep(0.05, length(Es.time)))

if (length(Es.time) > 0) {
    for (i in 1:length(Es.time)) {
        fit.paramnames <- c(fit.paramnames, paste("E", i, sep=""))
    }
}
