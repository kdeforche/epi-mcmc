library("deSolve")
system("R CMD SHLIB model-tv-cmp.cpp")
dyn.load("model-tv-cmp.so")

InvalidDataOffset <- 10000
Initial <- 1

G <- 5.2
Tinc <- 3

died_rate <- 0.007

calcNormalProfile <- function(mean, sd)
{
    kbegin = max(0, ceiling(mean - sd * 3))
    kend = max(kbegin + 1, floor(mean + sd * 3))

    result = NULL
    result$kbegin = -kend
    result$kend = -kbegin
    result$values = numeric(result$kend - result$kbegin)
    i = 1
    for (k in kbegin:kend) {
        result$values[i] = pnorm(k-0.5, mean, sd) -
            pnorm(k+0.5, mean, sd)
        i = i + 1
    }

    result$values = result$values / sum(result$values)
    result$values = rev(result$values)

    result
}

calcLogNormalProfile <- function(mean, sd)
{
    logm = log(mean) - 0.5*log((sd/mean)^2 + 1)
    logsd = sqrt(log((sd/mean)^2 + 1))

    kbegin = max(0, ceiling(mean - sd * 3))
    kend = max(kbegin + 1, floor(mean + sd * 3))

    result = NULL
    result$kbegin = -kend
    result$kend = -kbegin
    result$values = numeric(result$kend - result$kbegin)
    i = 1
    for (k in kbegin:kend) {
        result$values[i] = plnorm(k-0.5, logm, logsd) -
            plnorm(k+0.5, logm, sd=logsd)
        i = i + 1
    }

    result$values = result$values / sum(result$values)
    result$values = rev(result$values)

    result
}

calcGammaProfile <- function(mean, sd)
{
    shape = mean^2 / sd^2
    scale = sd^2 / mean
    
    kbegin = max(0, ceiling(mean - sd * 3))
    kend = max(kbegin + 1, floor(mean + sd * 3))

    result = NULL
    result$kbegin = -kend
    result$kend = -kbegin
    result$values = numeric(result$kend - result$kbegin)
    i = 1
    for (k in kbegin:kend) {
        result$values[i] = pgamma(k-0.5, shape=shape, scale=scale) -
            pgamma(k+0.5, shape=shape, scale=scale)
        i = i + 1
    }

    result$values = result$values / sum(result$values)
    result$values = rev(result$values)

    result
}

hospProfile <- calcGammaProfile
diedProfile <- calcGammaProfile

convolute <- function(values, i1, i2, profile)
{
    filter(values, profile$values, method="convolution", sides=1)[(i1 + profile$kend):(i2 + profile$kend)]
}

calculateModel <- function(params, period)
{
    Rt0 <- params[1]
    Rt1 <- params[2]
    Rt2 <- params[3]
    hosp_rate <- params[4]
    hosp_latency <- params[5]
    died_latency <- params[6]
    phs_morts <- params[7]
    phs <- params[8]
    HLsd <- params[9]
    DLsd <- params[10]

    Tinf <- (G - Tinc) * 2
    Tinft0 = Tinft1 = Tinft2 = Tinf
    
    if (phs > 0) { ## no physical distancing
        Rt01 = (Rt0 + Rt1) / 2
        Tinft01 = (Tinft0 + Tinft1) / 2

        Rt0 = Rt1 = Rt01
        Tinft0 = Tinft1 = Tinft01
    } 

    betat0 = Rt0 / Tinft0
    betat1 = Rt1 / Tinft1
    betat2 = Rt2 / Tinft2

    a <- 1 / Tinc

    ## https://www.imperial.ac.uk/media/imperial-college/medicine/mrc-gida/2020-06-08-COVID19-Report-26.pdf p.4 : gamma mean = 18.8, sd = 8.46
    hosp_cv_profile = hospProfile(hosp_latency, HLsd)
    died_cv_profile = diedProfile(died_latency, DLsd)

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
               phts = 1E10, ldts = 1E10, ldte = 1E10,
               betat0 = betat0, Tinft0 = Tinft0,
               betat1 = betat1, Tinft1 = Tinft1,
               betat2 = betat2, Tinft2 = Tinft2)

    Y <- c(S = N - Initial, E = Initial, I = 0, R = 0)

    times <- (padding + 1):(padding + period)

    out <- ode(Y, times, func = "derivs", parms = parms,
               jacfunc = "jac", dllname = "model-tv-cmp",
               initfunc = "initmod", nout = 2, outnames = c("Re", "Rt"))

    state$S[(padding + 1):(padding + period)] = out[,2]

    s2 <- convolute(state$S, padding + 1, padding + period, died_cv_profile)
    state$died[(padding + 1):(padding + period)] = (N - s2) * died_rate
    
    data_offset = InvalidDataOffset

    lds <- which(state$died > phs_morts)
    if (length(lds) > 0) {
        data_offset <- lds[1] - lockdown_offset

        parms <- c(N = N, a = a,
                   phts = data_offset + lockdown_offset + phs,
                   ldts = data_offset + lockdown_offset + max(phs, 0),
                   ldte = data_offset + lockdown_offset + lockdown_transition_period,
                   betat0 = betat0, Tinft0 = Tinft0,
                   betat1 = betat1, Tinft1 = Tinft1,
                   betat2 = betat2, Tinft2 = Tinft2)

        out <- ode(Y, times, func = "derivs", parms = parms,
                   jacfunc = "jac", dllname = "model-tv-cmp",
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
    params
}

invTransformParams <- function(posterior)
{
    posterior$Tinc = Tinc

    Tinf <- (G - Tinc) * 2
    posterior$Tinft0 = posterior$Tinft1 = posterior$Tinft2 = Tinf

    posterior$betat0 = posterior$Rt0 / posterior$Tinft0
    posterior$betat1 = posterior$Rt1 / posterior$Tinft1
    posterior$betat2 = posterior$Rt2 / posterior$Tinft2

    posterior$frRt0_1 = posterior$Rt1 / posterior$Rt0
    posterior$frRt1_2 = posterior$Rt2 / posterior$Rt1
    
    posterior
}

calclogp <- function(params) {
    Rt0 <- params[1]
    Rt1 <- params[2]
    Rt2 <- params[3]
    hosp_rate <- params[4]
    hosp_latency <- params[5]
    died_latency <- params[6]
    phs_morts <- params[7]
    phs <- params[8]
    HLsd <- params[9]
    DLsd <- params[10]

    logPriorP <- 0

    logPriorP <- logPriorP + dnorm(phs, mean=0, sd=10, log=T)
    logPriorP <- logPriorP + dnorm(Rt0 - Rt1, mean=0, sd=1, log=T)
    logPriorP <- logPriorP + dnorm(Rt1 - Rt2, mean=0, sd=1, log=T)
    logPriorP <- logPriorP + dnorm(HLsd, mean=5, sd=1, log=T)
    logPriorP <- logPriorP + dnorm(DLsd, mean=5, sd=1, log=T)

    logPriorP
}

## log likelihood function for fitting this model to observed data:
##   dhospi and dmorti
calclogl <- function(params, x) {
    state <<- calculateModel(params, FitTotalPeriod)

    if (state$offset == InvalidDataOffset)
        state$offset = 1

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

    result <- loglH + loglD
    
    if (it %% 1000 == 0) {
        print(params)
	print(c(it, result))
	graphs()
    }

    result
}

fit.paramnames <- c("Rt0", "Rt1", "Rt2", "HR", "HL", "DL", "phs_morts", "phs", "HLsd", "DLsd")
keyparamnames <- c("Rt0", "Rt1", "Rt2", "phs")
fitkeyparamnames <- c("Rt0", "Rt1", "Rt2", "phs")
init <- c(2.9, 0.9, 0.9, 0.02, 10, 21, total_deaths_at_lockdown, 0, 5, 5)

df_params <- data.frame(name = fit.paramnames,
                        min = c(0.1, 0.1, 0.1, 0.001, 5, 10, 0, -30, 2, 2),
                        max = c(8, 8, 8, 1, 25, 50,
                                max(dmort[length(dmort)] / 10, total_deaths_at_lockdown * 10),
                                30, 9, 9),
                        init = init)
