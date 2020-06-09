library("deSolve")
system("R CMD SHLIB model-vm-ode-cmp.cpp")
dyn.load("model-vm-ode-cmp.so")

InvalidDataOffset <- 10000
Initial <- 1

Tinf <- 8
Tinc <- 3
died_rate <- 0.007

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
    G0 <- params[3]
    Gt <- params[4]
    Tin0 <- params[5]
    Tint <- params[6]
    hosp_rate <- params[7]
    hosp_latency <- params[8]
    died_latency <- params[9]
    mort_lockdown_threshold <- params[10]
    lockdown_transition_start_var <- params[11]
    lockdown_transition_end_var <- params[12]

    Tis0 = Tinf
    K0 = Tin0 + 0.5 * Tinf
    Tina0 = (1 / (1 / Tin0 + 1 / Tinf))
    L0 = 0.5 * Tina0
    Rin0 = (K0 - (G0 - Tinc)) * R0 / (K0 - L0)
    Ris0 = R0 - Rin0
    betaIn0 = Rin0 / Tina0
    betaIs0 = Ris0 / Tis0
    beta0 = R0 / Tinf

    Tist = Tinf
    Kt = Tint + 0.5 * Tinf
    Tinat = (1 / (1 / Tint + 1 / Tinf))
    Lt = 0.5 * Tinat
    Rint = (Kt - (Gt - Tinc)) * Rt / (Kt - Lt)
    Rist = Rt - Rint
    betaInt = Rint / Tinat
    betaIst = Rist / Tist
    betat = Rt / Tinf
    
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
    state$Sr <- c()
    state$i <- padding + 1

    parms <- c(N = N, a = a, gamma = 1/Tinf,
               ldts = 1E10, ldte = 1E10,
               betaIn0 = betaIn0, betaIs0 = betaIs0, Tin0 = Tin0,
               betaInt = betaInt, betaIst = betaIst, Tint = Tint)

    Y <- c(S = N - Initial, E = Initial, In = 0, Is = 0, R = 0)

    times <- (padding + 1):(padding + period)

    out <- ode(Y, times, func = "derivs", parms = parms,
               jacfunc = "jac", dllname = "model-vm-ode-cmp",
               initfunc = "initmod", nout = 2, outnames = c("Re", "Rt"))

    state$S[(padding + 1):(padding + period)] = out[,2]

    s2 <- convolute(state$S, padding + 1, padding + period, died_cv_profile)
    state$died[(padding + 1):(padding + period)] = (N - s2) * died_rate
    
    data_offset = InvalidDataOffset          
    
    lds <- which(state$died > mort_lockdown_threshold)
    if (length(lds) > 0) {
        data_offset <- lds[1] - lockdown_offset - lockdown_transition_start_var
        parms <- c(N = N, a = a, gamma = 1/Tinf,
                   ldts = data_offset + lockdown_offset + lockdown_transition_start_var,
                   ldte = data_offset + lockdown_offset + lockdown_transition_period + lockdown_transition_end_var,
                   betaIn0 = betaIn0, betaIs0 = betaIs0, Tin0 = Tin0,
                   betaInt = betaInt, betaIst = betaIst, Tint = Tint)

        out <- ode(Y, times, func = "derivs", parms = parms,
               jacfunc = "jac", dllname = "model-vm-ode-cmp",
               initfunc = "initmod", nout = 2, outnames = c("Re", "Rt"))
    }

    state$S[(padding + 1):(padding + period)] = out[,2]
    state$E[(padding + 1):(padding + period)] = out[,3]
    state$In[(padding + 1):(padding + period)] = out[,4]
    state$Is[(padding + 1):(padding + period)] = out[,5]
    state$R[(padding + 1):(padding + period)] = out[,6]
    state$Re[(padding + 1):(padding + period)] = out[,7]
    state$Rt[(padding + 1):(padding + period)] = out[,8]

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
    result[7] = exp(params[7])

    result
}

invTransformParams <- function(posterior)
{
    posterior$Tinf = Tinf
    posterior$Tinc = Tinc

    posterior$Tis0 = posterior$Tinf
    K0 = posterior$Tin0 + 0.5 * posterior$Tinf
    Tina0 = (1 / (1 / posterior$Tin0 + 1 / posterior$Tinf))
    L0 = 0.5 * Tina0
    posterior$Rin0 = (K0 - (posterior$G0 - posterior$Tinc)) * posterior$R0 / (K0 - L0)
    posterior$Ris0 = posterior$R0 - posterior$Rin0
    posterior$betaIn0 = posterior$Rin0 / Tina0
    posterior$betaIs0 = posterior$Ris0 / posterior$Tis0

    posterior$Tinf.eff0 = 2 * (posterior$G0 - posterior$Tinc)
    posterior$beta.eff0 = posterior$R0 / posterior$Tinf.eff0
    
    posterior$Tist = posterior$Tinf
    Kt = posterior$Tint + 0.5 * posterior$Tinf
    Tinat = (1 / (1 / posterior$Tint + 1 / posterior$Tinf))
    Lt = 0.5 * Tinat
    posterior$Rint = (Kt - (posterior$Gt - posterior$Tinc)) * posterior$Rt / (Kt - Lt)
    posterior$Rist = posterior$Rt - posterior$Rint
    posterior$betaInt = posterior$Rint / Tinat
    posterior$betaIst = posterior$Rist / posterior$Tist

    posterior$Tinf.efft = 2 * (posterior$Gt - posterior$Tinc)
    posterior$beta.efft = posterior$Rt / posterior$Tinf.efft
    
    posterior$HR = exp(posterior$logHR)

    posterior
}

calclogp <- function(params) {
    R0 <- params[1]
    Rt <- params[2]
    G0 <- params[3]
    Gt <- params[4]
    Tin0 <- params[5]
    Tint <- params[6]
    hosp_rate <- params[7]
    hosp_latency <- params[8]
    died_latency <- params[9]
    mort_lockdown_threshold <- params[10]
    lockdown_transition_start_var <- params[11]
    lockdown_transition_end_var <- params[12]

    Tis0 = Tinf
    K0 = Tin0 + 0.5 * Tinf
    Tina0 = (1 / (1 / Tin0 + 1 / Tinf))
    L0 = 0.5 * Tina0
    Rin0 = (K0 - (G0 - Tinc)) * R0 / (K0 - L0)
    Ris0 = R0 - Rin0
    betaIn0 = Rin0 / Tina0
    betaIs0 = Ris0 / Tis0
    beta0 = R0 / Tinf

    Tist = Tinf
    Kt = Tint + 0.5 * Tinf
    Tinat = (1 / (1 / Tint + 1 / Tinf))
    Lt = 0.5 * Tinat
    Rint = (Kt - (Gt - Tinc)) * Rt / (Kt - Lt)
    Rist = Rt - Rint
    betaInt = Rint / Tinat
    betaIst = Rist / Tist
    betat = Rt / Tinf

    if (betaIn0 < 0 || betaInt < 0 || betaIs0 < 0 || betaIst < 0) {
        return(-Inf)
    }
    
    if (betaIst > betaInt) {
        return(-Inf)
    }

    if (betaIs0 > betaIn0) {
        return(-Inf)
    }

    logPriorP <- 0

    logPriorP <- logPriorP + dnorm(lockdown_transition_start_var, mean=0, sd=3, log=T)

    if (lockdown_transition_end_var > lockdown_transition_start_var)
        return(-Inf)
    
    logPriorP <- logPriorP + dnorm(lockdown_transition_end_var, mean=0, sd=3, log=T)

    logPriorP <- logPriorP + dnorm(died_latency, mean=21, sd=4, log=T)
    logPriorP <- logPriorP + dnorm(G0, mean=5, sd=1, log=T)
    logPriorP <- logPriorP + dnorm(G0 - Gt, mean=0, sd=1.5, log=T)

    logPriorP
}

calclogl <- function(params) {
    R0 <- params[1]
    Rt <- params[2]
    G0 <- params[3]
    Gt <- params[4]
    Tin0 <- params[5]
    Tint <- params[6]
    hosp_rate <- params[7]
    hosp_latency <- params[8]
    died_latency <- params[9]
    mort_lockdown_threshold <- params[10]
    lockdown_transition_start_var <- params[11]
    lockdown_transition_end_var <- params[12]

    state <<- calculateModel(params, FitTotalPeriod)

    if (state$offset == InvalidDataOffset)
        state$offset = 1

    total_deaths_at_lockdown <-
        dmort[max(1, lockdown_offset + round(lockdown_transition_start_var))]

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
	print(c(it, result))
	graphs()
    }

    result
}

fit.paramnames <- c("R0", "Rt", "G0", "Gt", "Tin0", "Tint",
                    "logHR", "HL", "DL", "lockdownmort", "ldtsv", "ldtev")
keyparamnames <- c("betaIn0", "betaInt", "betaIs0", "betaIst", "R0", "Rt", "G0", "Gt",
                   "Rin0", "Ris0")
fitkeyparamnames <- c("R0", "Rt", "G0", "Gt", "Tin0", "Tint")
init <- c(2.9, 0.9, 5, 5, 5, 5, log(0.02), 10, 20, total_deaths_at_lockdown, 0, 0)

df_params <- data.frame(name = fit.paramnames,
                        min = c(0.1, 0.1, Tinc + 0.2, Tinc + 0.2, 0.2, 0.2, log(0.001), 5, 5, 0,
                                -lockdown_offset, -lockdown_offset),
                        max = c(8, 8, 8, 8, 7.8, 7.8, log(0.5), 30, 40,
                                max(dmort[length(dmort)] / 10, total_deaths_at_lockdown * 10),
                                40, 40),
                        init = init)
