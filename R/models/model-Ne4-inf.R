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

## Kleinere tijdseenheid: 0.5d, iterate here more than once if necessary:
## 0.5, 0.33, 0.25, ...

simstep <- function(state, N, betaIn, betaIs, a, Tin, Tis, k)
{
    i = state$i

    gamma1 = 1/Tin
    gamma2 = 1/Tis

    loops = ceiling(gamma1)
    Ts = 1/loops

    S = state$S[i]
    E = state$E[i]
    In = state$In[i]
    Is = state$Is[i]
    R = state$R[i]

    day_infected <- 0

    for (l in 1:loops) {
        EIR = N - S
        Ne = max(EIR, k * N)

        got_infected = (betaIn * In + betaIs * Is) / N * (Ne - EIR)
        got_infectious = a * E
        got_isolated = gamma1 * In
        got_removed = gamma2 * Is

        ##print(c(i, Tin, Tis, got_infected, got_infectious, got_isolated, got_removed))
    
        deltaS = -got_infected
        deltaE = got_infected - got_infectious
        deltaIn = got_infectious - got_isolated
        deltaIs = got_isolated - got_removed
        deltaR = got_removed

        S = S + deltaS * Ts
        E = E + deltaE * Ts
        In = In + deltaIn * Ts
        Is = Is + deltaIs * Ts
        R = R + deltaR * Ts

        if (l == 1)
            day_infected <- got_infected
    }

    ## only based on first 'loop iteration'
    state$Re[i + 1] = day_infected / ((state$In[i] + state$Is[i]) / (Tin + Tis))
    state$Rt[i + 1] = state$Re[i + 1] * N / (Ne - EIR)
    state$S[i + 1] = S
    state$E[i + 1] = E
    state$In[i + 1] = In
    state$Is[i + 1] = Is
    state$R[i + 1] = R

    i = i + 1
    state$i = i

    state
}

##
## Parameters are a function of time:
##  t < lockdown_offset : par0
##  t > lockdown_offset + lockdown_transition_period : part
##  t > Es.time[i] : Es[i] * part + (1 - Es[i]) * par0
##

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
    betaIn <- params[1]
    betaIs <- params[2]
    hosp_rate <- params[3]
    died_rate <- params[4]
    hosp_latency <- params[5]
    died_latency <- params[6]
    Tinf <- params[7]
    Tinc <- params[8]
    mort_lockdown_threshold <- params[9]
    k0 <- params[10]
    kt <- params[11]
    Tin0 <- params[12]
    Tint <- params[13]
    Es <- tail(params, n=-13)

    h0 <- 1
    a <- 1 / Tinc

    hosp_cv_profile = calcConvolveProfile(-hosp_latency, 5)
    died_cv_profile = calcConvolveProfile(-died_latency, 5)

    padding = max(-hosp_cv_profile$kbegin, -died_cv_profile$kbegin) + 1

    state <- NULL
    state$Re <- rep(0, padding + period)
    state$Rt <- rep(0, padding + period)
    state$S <- rep(N-Initial, padding + period)
    state$E <- rep(Initial, padding + period)
    state$In <- rep(0, padding + period)
    state$Is <- rep(0, padding + period)
    state$R <- rep(0, padding + period)
    state$hosp <- rep(0, padding + period)
    state$died <- rep(0, padding + period)
    state$i <- padding + 1

    data_offset = InvalidDataOffset

    if (Tin0 > 0 && Tint > 0 && Tinc > 1) {
        for (i in (padding + 1):(padding + period)) {
            k = calcpar(i - data_offset, k0, kt, Es)
            Tin = calcpar(i - data_offset, Tin0, Tint, Es)
            Tis = Tinf - Tin
            
            state <- simstep(state, N, betaIn, betaIs, a, Tin, Tis, k)

            s = convolute(state$S, i, hosp_cv_profile)
            state$hosp[i] <- (N - s) * hosp_rate

            ## died:
            ##  from those that (ever) became infectious           
            r = convolute(state$In + state$Is + state$R, i, died_cv_profile)
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

    result
}

invTransformParams <- function(posterior)
{
    posterior$IFR = exp(posterior$logHR + posterior$logHRDR)
    posterior$HR = exp(posterior$logHR)

    ## Additional quantitites of interest
    ## average duration one stays infecting * infection rate
    ## + average duration one stays isolated * isolation infection rate
    posterior$R0 = posterior$betaIn * posterior$Tin0 + posterior$betaIs * (posterior$Tinf - posterior$Tin0)
    posterior$Rt = posterior$betaIn * posterior$Tint + posterior$betaIs * (posterior$Tinf - posterior$Tint)

    posterior
}

## log likelihood function for fitting this model to observed data:
##   dhospi and dmorti
calclogl <- function(params) {
    betaIn <- params[1]
    betaIs <- params[2]
    hosp_rate <- params[3]
    died_rate <- params[4]
    hosp_latency <- params[5]
    died_latency <- params[6]
    Tinf <- params[7]
    Tinc <- params[8]
    mort_lockdown_threshold <- params[9]
    k0 <- params[10]
    kt <- params[11]
    Tin0 <- params[12]
    Tint <- params[13]
    Es <- tail(params, n=-13)

    if (betaIn < 1E-10 || betaIs < 1E-10) {
        return(-Inf)
    }

    if (Tin0 < 0 || Tin0 > Tinf - 1) {
        return(-Inf)
    }

    if (Tint < 0 || Tint > Tin0) {
        return(-Inf)
    }

    if (k0 <= 0 || k0 > 1) {
        return(-Inf)
    }

    if (kt <= 0 || kt > 1) {
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

    ## logPriorP <- logPriorP + dnorm(beta, 1, 0.22, log=T)
    logPriorP <- logPriorP + dbeta(k0, 9, 1, log=T)
    logPriorP <- logPriorP + dbeta(kt, 4, 1, log=T)
    logPriorP <- logPriorP + dnorm(betaIs, 0, 0.5, log=T)

    ## estBetaParams(0.0066, 0.005^2)
    ## logPriorP <- logPriorP + dbeta(died_rate, 1.7243, 259.53, log=T)
    ## estBetaParams(0.0066, 0.002^2)
    ##logPriorP <- logPriorP + dbeta(died_rate, 10.811, 1627.298, log=T)
    ## estBetaParams(0.0066, 0.001^2)
    logPriorP <- logPriorP + dbeta(died_rate, 43.2659, 6512.174, log=T)

    logPriorP <- logPriorP + dnorm(hosp_latency, mean=10, sd=20, log=T)
    logPriorP <- logPriorP + dnorm(died_latency, mean=10, sd=20, log=T)

    logPriorP <- logPriorP + dnorm(Tinc, mean=5, sd=1, log=T)
    logPriorP <- logPriorP + dnorm(Tinf, mean=8, sd=1, log=T)

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

    if (dend < 1 || dstart < 1) {
        print(c(dstart, dend))
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

fit.paramnames <- c("betaIn", "betaIs", "logHR", "logHRDR", "HL", "DL",
                    "Tinf", "Tinc", "lockdownmort", "k0", "kt", "Tin0", "Tint")
keyparamnames <- c("betaIn", "betaIs", "IFR", "Tinf", "Tinc", "R0", "Rt", "k0", "kt")
fitkeyparamnames <- c("betaIn", "betaIs", "logHRDR", "Tinf", "Tinc", "k0", "kt", "Tin0", "Tint")

init <- c(1, 1, log(0.02), log(0.3), 10, 20, 8, 5, total_deaths_at_lockdown, 0.9, 0.5, 7, 2,
          rep(0.9, length(Es.time)))
scales <- c(0.3, 0.3, 0.05, 0.1, 1, 1, 1, 1, total_deaths_at_lockdown / 20, 0.3, 0.3, 1, 1,
            rep(0.05, length(Es.time)))

if (length(Es.time) > 0) {
    for (i in 1:length(Es.time)) {
        fit.paramnames <- c(fit.paramnames, paste("E", i, sep=""))
    }
}
