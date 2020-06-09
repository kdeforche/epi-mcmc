library("deSolve")
system("R CMD SHLIB model-simple-ode-cmp.cpp")
dyn.load("model-simple-ode-cmp.so")

InvalidDataOffset <- 10000
Initial <- 1

Tinc <- 3
died_rate <- 0.007

if (!exists("Es.time")) {
    Es.time <- c()
}

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

convolute <- function(values, i1, i2, profile)
{
    filter(values, profile$values, method="convolution", sides=1)[(i1 + profile$kend):(i2 + profile$kend)]
}

calculateModel <- function(params, period)
{
    R0 <- params[1]
    Rt <- params[2]
    Tinf0 <- params[3]
    Tinft <- params[4]
    hosp_rate <- params[5]
    hosp_latency <- params[6]
    died_latency <- params[7]
    mort_lockdown_threshold <- params[8]
    Es <- tail(params, n=-8)

    beta0 = R0 / Tinf0
    betat = Rt / Tinft
    
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
    state$Sr <- c()
    state$i <- padding + 1

    parms <- c(N = N, a = a,
               ldts = 1E10, ldte = 1E10,
               beta0 = beta0, Tinf0 = Tinf0,
               betat = betat, Tinft = Tinft)

    Y <- c(S = N - Initial, E = Initial, I = 0, R = 0)

    times <- (padding + 1):(padding + period)

    out <- ode(Y, times, func = "derivs", parms = parms,
               jacfunc = "jac", dllname = "model-simple-ode-cmp",
               initfunc = "initmod", nout = 2, outnames = c("Re", "Rt"))

    state$S[(padding + 1):(padding + period)] = out[,2]

    s2 <- convolute(state$S, padding + 1, padding + period, died_cv_profile)
    state$died[(padding + 1):(padding + period)] = (N - s2) * died_rate
    
    data_offset = InvalidDataOffset

    lds <- which(state$died > mort_lockdown_threshold)
    if (length(lds) > 0) {
        data_offset <- lds[1] - lockdown_offset
        parms <- c(N = N, a = a,
                   ldts = data_offset + lockdown_offset,
                   ldte = data_offset + lockdown_offset + lockdown_transition_period,
                   beta0 = beta0, Tinf0 = Tinf0,
                   betat = betat, Tinft = Tinft)

        out <- ode(Y, times, func = "derivs", parms = parms,
                   jacfunc = "jac", dllname = "model-simple-ode-cmp",
                   initfunc = "initmod", nout = 2, outnames = c("Re", "Rt"))
    }

    state$S[(padding + 1):(padding + period)] = out[,2]
    state$E[(padding + 1):(padding + period)] = out[,3]
    state$I[(padding + 1):(padding + period)] = out[,4]
    state$R[(padding + 1):(padding + period)] = out[,5]
    state$Re[(padding + 1):(padding + period)] = out[,6]
    state$Rt[(padding + 1):(padding + period)] = out[,7]

    s1 <- convolute(state$S, padding + 1, padding + period, hosp_cv_profile)
    state$hosp[(padding + 1):(padding + period)] = (N - s1) * hosp_rate
    
    s2 <- convolute(state$S, padding + 1, padding + period, died_cv_profile)
    state$died[(padding + 1):(padding + period)] = (N - s2) * died_rate

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
    result[5] = exp(params[5])

    result
}

invTransformParams <- function(posterior)
{
    posterior$Tinc = Tinc

    posterior$beta0 = posterior$R0 / posterior$Tinf0
    posterior$betat = posterior$Rt / posterior$Tinft

    posterior$HR = exp(posterior$logHR)

    posterior
}

calclogp <- function(params) {
    R0 <- params[1]
    Rt <- params[2]
    Tinf0 <- params[3]
    Tinft <- params[4]
    hosp_rate <- params[5]
    hosp_latency <- params[6]
    died_latency <- params[7]
    mort_lockdown_threshold <- params[8]
    Es <- tail(params, n=-8)

    if (R0 < 0 || Rt < 0) {
        return(-Inf)
    }
    
    if (Tinf0 < 0.2 || Tinf0 > 13.8) {
        return(-Inf)
    }

    if (Tinft < 0.2 || Tinft > 13.8) {
        return(-Inf)
    }

    beta0 = R0 / Tinf0
    G0 = Tinc + Tinf0/2
    betat = Rt / Tinft
    Gt = Tinc + Tinft/2

    if (beta0 < 0 || betat < 0) {
        return(-Inf)
    }
    
    if (hosp_latency < 0 || hosp_latency > 30) {
        return(-Inf)
    }

    if (died_latency < 0 || died_latency > 30) {
        return(-Inf)
    }

    logPriorP <- 0

    ##logPriorP <- logPriorP + dnorm(R0, mean=2, sd=5, log=T)
    ##logPriorP <- logPriorP + dnorm(Rt, mean=2, sd=5, log=T)
    ##logPriorP <- logPriorP + dnorm(Tinf0, mean=4, sd=5, log=T)
    ##logPriorP <- logPriorP + dnorm(Tinf0, mean=4, sd=5, log=T)
    
    logPriorP <- logPriorP + dnorm(hosp_latency, mean=10, sd=20, log=T)
    logPriorP <- logPriorP + dnorm(died_latency, mean=21, sd=4, log=T)

    logPriorP <- logPriorP + dnorm(G0, mean=5, sd=1, log=T)
    logPriorP <- logPriorP + dnorm(G0 - Gt, mean=0, sd=1.5, log=T)

    logPriorP
}

## log likelihood function for fitting this model to observed data:
##   dhospi and dmorti
calclogl <- function(params, x) {
    R0 <- params[1]
    Rt <- params[2]
    Tinf0 <- params[3]
    Tinft <- params[4]
    hosp_rate <- params[5]
    hosp_latency <- params[6]
    died_latency <- params[7]
    mort_lockdown_threshold <- params[8]

    beta0 = R0 / Tinf0
    G0 = Tinc + Tinf0/2
    betat = Rt / Tinft
    Gt = Tinc + Tinft/2

    state <<- calculateModel(params, FitTotalPeriod)

    if (state$offset == InvalidDataOffset)
        state$offset = 1

    loglLD <- dnbinom(total_deaths_at_lockdown, mu=pmax(0.1, mort_lockdown_threshold),
                      size=mort_nbinom_size, log=T)

    dstart <- state$offset
    dend <- state$offset + length(dhospi) - 1

    if (dend > length(state$hospi)) {
        ##print("=========================== Increase FitTotalPeriod ===================")
        dend <- length(state$hospi)
        dstart <- dend - length(dhospi) + 1
    }

    di <- dhospi

    if (dstart < 1) {
        ##print("=========================== Increase Padding ? ===================")
        dstart <- 1
        dend <- dstart + length(dhospi) - 1
    }
    
    loglH <- sum(dnbinom(dhospi,
                         mu=pmax(0.1, state$hospi[dstart:dend]),
                         size=hosp_nbinom_size, log=T))

    dstart <- state$offset
    dend <- state$offset + length(dmorti) - 1

    if (dend > length(state$deadi)) {
        ##print("=========================== Increase FitTotalPeriod ===================")
        dend <- length(state$deadi)
        dstart <- dend - length(dmorti) + 1
    }

    if (dstart < 1) {
        ##print("=========================== Increase Padding ? ===================")
        dstart <- 1
        dend <- dstart + length(dmorti) - 1
    }

    loglD <- sum(dnbinom(dmorti,
                         mu=pmax(0.1, state$deadi[dstart:dend]), size=mort_nbinom_size, log=T))

    it <<- it + 1

    result <- loglLD + loglH + loglD
    
    if (it %% 1000 == 0) {
        print(params)
        print(c(R0, Rt, Tinf0, Tinft, beta0, betat))
	print(c(it, result))
	graphs()
    }

    result
}

fit.paramnames <- c("R0", "Rt", "Tinf0", "Tinft",
                    "logHR", "HL", "DL", "lockdownmort")
keyparamnames <- c("beta0", "betat", "R0", "Rt", "Tinf0", "Tinft")
fitkeyparamnames <- c("R0", "Rt", "Tinf0", "Tinft")
init <- c(2.9, 0.9, 3, 3, log(0.02), 10, 20, total_deaths_at_lockdown)
scales <- c(1, 1, 1, 1, 0.05, 1, 1, total_deaths_at_lockdown / 20)

df_params <- data.frame(name = fit.paramnames,
                        min = c(0.1, 0.1, 0.1, 0.1, log(0.001), 5, 5, 0),
                        max = c(8, 2, 8, 8, log(0.5), 30, 30,
                                max(dmort[length(dmort)] / 10, total_deaths_at_lockdown * 10)),
                        init = init)
