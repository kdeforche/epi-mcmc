library("deSolve")
system("R CMD SHLIB model-agef.cpp")
dyn.load("model-agef.so")

InvalidDataOffset <- 10000
Initial <- 1

G <- 4.7
Tinc1 <- 3
Tinc2 <- 1
y.died_rate <- y.ifr[1]
o.died_rate <- o.ifr[1]
a1 <- 1/Tinc1
a2 <- 1/Tinc2
gamma <- 1/((G - (Tinc1 + Tinc2)) * 2)

calcGammaProfile <- function(mean, sd)
{
    scale = sd^2 / mean
    shape = mean / scale
    
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

convolute <- function(values, i1, i2, profile)
{
    filter(values, profile$values, method="convolution", sides=1)[(i1 + profile$kend):(i2 + profile$kend)]
}

caseProfile <- calcGammaProfile
diedProfile <- calcGammaProfile

calculateModel <- function(params, period)
{
    betay0 <- params[1]
    betao0 <- params[2]
    betayo0 <- params[3]
    betay1 <- params[4]
    betao1 <- params[5]
    betayo1 <- params[6]
    betay2 <- params[7]
    betao2 <- params[8]
    betayo2 <- params[9]
    ycase_rate <- params[10]
    ycase_latency <- params[11]
    yhosp_rate <- params[12]
    yhosp_latency <- params[13]
    ydied_latency <- params[14]
    ocase_rate <- params[15]
    ocase_latency <- params[16]
    ohosp_rate <- params[17]
    ohosp_latency <- params[18]
    odied_latency <- params[19]
    t0_morts <- params[20]
    t0o <- params[21]
    CLsd <- params[22]                  
    HLsd <- params[23]
    DLsd <- params[24]
    t3o <- params[25] ## start of increase, 6 June
    betay3 <- betay2
    betao3 <- betao2
    betayo3 <- betayo2
    t4o <- params[26] ## end of increase, 6 Juli
    betay4 <- params[27]
    betao4 <- params[28]
    betayo4 <- params[29]
    betay5 <- betay4 ## new lockdown, 27 Juli
    betao5 <- betao4
    betayo5 <- betayo4
    t6o <- params[30]
    betay6 <- params[31] ## end of lockdown transition 27 Juli + t6o
    betao6 <- params[32]
    betayo6 <- params[33]
    betay7 <- params[34] * betay6 ## start of school Sept 1
    betayo7 <- params[35] * betayo6
    betao7 <- params[36] * betao6
    
    ## convolution profiles
    y.case_cv_profile = caseProfile(ycase_latency, CLsd)
    o.case_cv_profile = caseProfile(ocase_latency, CLsd)
    y.hosp_cv_profile = caseProfile(yhosp_latency, HLsd)
    o.hosp_cv_profile = caseProfile(ohosp_latency, HLsd)
    y.died_cv_profile = diedProfile(ydied_latency, DLsd)
    o.died_cv_profile = diedProfile(odied_latency, DLsd)

    padding = max(-y.case_cv_profile$kbegin, -y.died_cv_profile$kbegin,
                  -o.case_cv_profile$kbegin, -o.died_cv_profile$kbegin) + 1

    state <- NULL

    state$y.S <- rep(y.N - Initial, padding + period)
    state$y.E1 <- rep(Initial, padding + period)
    state$y.E2 <- rep(0, padding + period)
    state$y.I <- rep(0, padding + period)
    state$y.R <- rep(0, padding + period)
    state$y.casei <- rep(0, padding + period)
    state$y.hospi <- rep(0, padding + period)
    state$y.deadi <- rep(0, padding + period)

    state$o.S <- rep(o.N, padding + period)
    state$o.E1 <- rep(0, padding + period)
    state$o.E2 <- rep(0, padding + period)
    state$o.I <- rep(0, padding + period)
    state$o.R <- rep(0, padding + period)
    state$o.casei <- rep(0, padding + period)
    state$o.hospi <- rep(0, padding + period)
    state$o.deadi <- rep(0, padding + period)

    state$Re <- rep(0, padding + period)
    state$Rt <- rep(0, padding + period)
    state$y.Re <- rep(0, padding + period)
    state$o.Re <- rep(0, padding + period)

    state$i <- padding + 1
    
    parms <- c(Ny = y.N, No = o.N,
               a1 = a1, a2 = a2, gamma = gamma,
               t0 = 1E10, t1 = 1E10, t2 = 1E10,
               betay0 = betay0, betao0 = betao0, betayo0 = betayo0,
               betay1 = betay1, betao1 = betao1, betayo1 = betayo1,
               betay2 = betay2, betao2 = betao2, betayo2 = betayo2,
               t3 = 1E10, betay3 = betay3, betao3 = betao3, betayo3 = betayo3,
               t4 = 1E10, betay4 = betay4, betao4 = betao4, betayo4 = betayo4,
               t5 = 1E10, betay5 = betay5, betao5 = betao5, betayo5 = betayo5,
               t6 = 1E10, betay6 = betay6, betao6 = betao6, betayo6 = betayo6,
               t7 = 1E10, betay7 = betay7, betao7 = betao7, betayo7 = betayo7)

    Y <- c(Sy = y.N - Initial, E1y = Initial, E2y = 0, Iy = 0, Ry = 0,
           So = o.N, E1o = 0, E2o = 0, Io = 0, Ro = 0)

    times <- (padding + 1):(padding + period)

    out <- ode(Y, times, func = "derivs", parms = parms,
               dllname = "model-agef",
               initfunc = "initmod", nout = 4, outnames = c("Re", "Rt", "y.Re", "o.Re"))

    state$y.S[(padding + 1):(padding + period)] = out[,2]
    state$o.S[(padding + 1):(padding + period)] = out[,7]

    s2 <- convolute(state$y.S, padding + 1, padding + period, y.died_cv_profile)
    state$y.died[(padding + 1):(padding + period)] = (y.N - s2) * y.died_rate

    s2 <- convolute(state$o.S, padding + 1, padding + period, o.died_cv_profile)
    state$o.died[(padding + 1):(padding + period)] = (o.N - s2) * o.died_rate

    state$died = state$y.died + state$o.died

    data_offset = InvalidDataOffset

    lds <- which(state$died > t0_morts)
    if (length(lds) > 0) {
        data_offset <- lds[1] - lockdown_offset

        t2 <- data_offset + lockdown_offset + lockdown_transition_period
        t3 <- t2 + t3o
        t4 <- t3 + t4o

        t5 <- data_offset + d5        
        t6 <- t5 + t6o
        t7 <- data_offset + d7

        parms <- c(Ny = y.N, No = o.N,
                   a1 = a1, a2 = a2, gamma = gamma,
                   t0 = data_offset + lockdown_offset + t0o,
                   t1 = data_offset + lockdown_offset + max(t0o, 0),
                   t2 = t2,
                   betay0 = betay0, betao0 = betao0, betayo0 = betayo0,
                   betay1 = betay1, betao1 = betao1, betayo1 = betayo1,
                   betay2 = betay2, betao2 = betao2, betayo2 = betayo2,
                   t3 = t3, betay3 = betay3, betao3 = betao3, betayo3 = betayo3,
                   t4 = t4, betay4 = betay4, betao4 = betao4, betayo4 = betayo4,
                   t5 = t5, betay5 = betay5, betao5 = betao5, betayo5 = betayo5,
                   t6 = t6, betay6 = betay6, betao6 = betao6, betayo6 = betayo6,
                   t7 = t7, betay7 = betay7, betao7 = betao7, betayo7 = betayo7)

        out <- ode(Y, times, func = "derivs", parms = parms,
                   dllname = "model-agef",
                   initfunc = "initmod", nout = 4, outnames = c("Re", "Rt", "y.Re", "o.Re"))
    }

    state$y.S[(padding + 1):(padding + period)] = out[,2]
    state$y.E[(padding + 1):(padding + period)] = out[,3] + out[,4]
    state$y.I[(padding + 1):(padding + period)] = out[,5]
    state$y.R[(padding + 1):(padding + period)] = out[,6]
    state$o.S[(padding + 1):(padding + period)] = out[,7]
    state$o.E[(padding + 1):(padding + period)] = out[,8]
    state$o.I[(padding + 1):(padding + period)] = out[,9] + out[,10]
    state$o.R[(padding + 1):(padding + period)] = out[,11]
    state$Re[(padding + 1):(padding + period)] = out[,12]
    state$Rt[(padding + 1):(padding + period)] = out[,13]
    state$y.Re[(padding + 1):(padding + period)] = out[,14]
    state$o.Re[(padding + 1):(padding + period)] = out[,15]

    ## y.case
    s2 <- y.N - convolute(state$y.S, padding + 1, padding + period, y.case_cv_profile)
    s2i <- c(s2[1], diff(s2))
    if (data_offset != InvalidDataOffset) {
        ## IFR is based on death, IHR is shifted in time (DL - CL) earlier
        t1 <- data_offset + round(-ydied_latency + ycase_latency)
        t2 <- min(padding + period, t1 + length(y.ifr) - 1)
        if (t1 > padding && t2 > t1) {
            state$y.casei[(padding + 1):t1] = s2i[1:(t1 - padding)] * pmin(1, ycase_rate * y.ifr[1])
            state$y.casei[t1:t2] = s2i[(t1 - padding):(t2 - padding)] * pmin(1, ycase_rate * y.ifr[1:(t2 - t1 + 1)])
            state$y.casei[t2:(padding + period)] = s2i[(t2 - padding):period] * pmin(1, ycase_rate * y.ifr[length(y.ifr)])
        }
    } else {
        state$y.casei[(padding + 1):(padding + period)] = s2i * ycase_rate * y.died_rate
    }
    state$y.case = cumsum(state$y.casei)

    ## y.hosp
    s2 <- y.N - convolute(state$y.S, padding + 1, padding + period, y.hosp_cv_profile)
    s2i <- c(s2[1], diff(s2))
    if (data_offset != InvalidDataOffset) {
        ## IFR is based on death, IHR is shifted in time (DL - HL) earlier
        t1 <- data_offset + round(-ydied_latency + yhosp_latency)
        t2 <- min(padding + period, t1 + length(y.ifr) - 1)
        if (t1 > padding && t2 > t1) {
            state$y.hospi[(padding + 1):t1] = s2i[1:(t1 - padding)] * yhosp_rate * y.hr[1]
            state$y.hospi[t1:t2] = s2i[(t1 - padding):(t2 - padding)] * yhosp_rate * y.hr[1:(t2 - t1 + 1)]
            state$y.hospi[t2:(padding + period)] = s2i[(t2 - padding):period] * yhosp_rate * y.hr[length(y.hr)]
        }
    } else {
        state$y.hospi[(padding + 1):(padding + period)] = s2i * yhosp_rate * y.hr[1]
    }

    ## y.dead
    s2 <- y.N - convolute(state$y.S, padding + 1, padding + period, y.died_cv_profile)
    s2i <- c(s2[1], diff(s2))
    if (data_offset != InvalidDataOffset) {
        t1 <- data_offset
        t2 <- min(padding + period, t1 + length(y.ifr) - 1)
        if (t1 > padding && t2 > t1) {
            state$y.deadi[(padding + 1):t1] = s2i[1:(t1 - padding)] * y.ifr[1]
            state$y.deadi[t1:t2] = s2i[(t1 - padding):(t2 - padding)] * y.ifr[1:(t2 - t1 + 1)]
            state$y.deadi[t2:(padding + period)] = s2i[(t2 - padding):period] * y.ifr[length(y.ifr)]
        }
    } else {
        state$y.deadi[(padding + 1):(padding + period)] = s2i * y.died_rate
    }
    state$y.died = cumsum(state$y.deadi)

    ## o.case
    s2 <- o.N - convolute(state$o.S, padding + 1, padding + period, o.case_cv_profile)
    s2i <- c(s2[1], diff(s2))
    if (data_offset != InvalidDataOffset) {
        t1 <- data_offset + round(-ydied_latency + ocase_latency)
        t2 <- min(padding + period, t1 + length(o.ifr) - 1)
        if (t1 > padding && t2 > t1) {
            state$o.casei[(padding + 1):t1] = s2i[1:(t1 - padding)] * pmin(1, ocase_rate * o.ifr[1])
            state$o.casei[t1:t2] = s2i[(t1 - padding):(t2 - padding)] * pmin(1, ocase_rate * o.ifr[1:(t2 - t1 + 1)])
            state$o.casei[t2:(padding + period)] = s2i[(t2 - padding):period] * pmin(1, ocase_rate * o.ifr[length(o.ifr)])
        }
    } else {
        state$o.casei[(padding + 1):(padding + period)] = s2i * ocase_rate * o.died_rate
    }
    state$o.case = cumsum(state$o.casei)

    ## o.hosp
    s2 <- o.N - convolute(state$o.S, padding + 1, padding + period, o.hosp_cv_profile)
    s2i <- c(s2[1], diff(s2))
    if (data_offset != InvalidDataOffset) {
        t1 <- data_offset + round(-ydied_latency + ohosp_latency)
        t2 <- min(padding + period, t1 + length(o.ifr) - 1)
        if (t1 > padding && t2 > t1) {
            state$o.hospi[(padding + 1):t1] = s2i[1:(t1 - padding)] * ohosp_rate * o.hr[1]
            state$o.hospi[t1:t2] = s2i[(t1 - padding):(t2 - padding)] * ohosp_rate * o.hr[1:(t2 - t1 + 1)]
            state$o.hospi[t2:(padding + period)] = s2i[(t2 - padding):period] * ohosp_rate * o.hr[length(o.hr)]
        }
    } else {
        state$o.hospi[(padding + 1):(padding + period)] = s2i * ohosp_rate * o.hr[1]
    }

    ## o.dead
    s2 <- o.N - convolute(state$o.S, padding + 1, padding + period, o.died_cv_profile)
    s2i <- c(s2[1], diff(s2))
    if (data_offset != InvalidDataOffset) {
        t1 <- data_offset
        t2 <- min(padding + period, t1 + length(o.ifr) - 1)
        if (t1 > padding && t2 > t1) {
            state$o.deadi[(padding + 1):t1] = s2i[1:(t1 - padding)] * o.ifr[1]
            state$o.deadi[t1:t2] = s2i[(t1 - padding):(t2 - padding)] * o.ifr[1:(t2 - t1 + 1)]
            state$o.deadi[t2:(padding + period)] = s2i[(t2 - padding):period] * o.ifr[length(o.ifr)]
        }
    } else {
        state$o.deadi[(padding + 1):(padding + period)] = s2i * o.died_rate
    }
    state$o.died = cumsum(state$o.deadi)

    state$padding <- padding
    state$offset <- data_offset

    state
}

transformParams <- function(params)
{
    params
}

invTransformParams <- function(posterior)
{
    posterior$Tinf <- 1/gamma

    posterior$y.Rt0 = posterior$betay0 / gamma
    posterior$o.Rt0 = posterior$betao0 / gamma
    posterior$yo.Rt0 = posterior$betayo0 / gamma

    posterior$y.Rt1 = posterior$betay1 / gamma
    posterior$o.Rt1 = posterior$betao1 / gamma
    posterior$yo.Rt1 = posterior$betayo1 / gamma

    posterior$y.Rt2 = posterior$betay2 / gamma
    posterior$o.Rt2 = posterior$betao2 / gamma
    posterior$yo.Rt2 = posterior$betayo2 / gamma

    posterior$y.Rt4 = posterior$betay4 / gamma
    posterior$o.Rt4 = posterior$betao4 / gamma
    posterior$yo.Rt4 = posterior$betayo4 / gamma

    posterior$y.Rt6 = posterior$betay6 / gamma
    posterior$o.Rt6 = posterior$betao6 / gamma
    posterior$yo.Rt6 = posterior$betayo6 / gamma

    posterior$betay7 = posterior$betay6 * posterior$fy7
    posterior$betao7 = posterior$betao6 * posterior$fo7
    posterior$betayo7 = posterior$betayo6 * posterior$fyo7

    posterior$y.Rt7 = posterior$betay7 / gamma
    posterior$o.Rt7 = posterior$betao7 / gamma
    posterior$yo.Rt7 = posterior$betayo7 / gamma

    posterior
}

calcNominalState <- function(state)
{
    state$case <- state$y.case + state$o.case
    state$died <- state$y.died + state$o.died
    state$casei <- state$y.casei + state$o.casei
    state$deadi <- state$y.deadi + state$o.deadi
    state$hospi <- state$y.hospi + state$o.hospi

    state
}

## log likelihood function for fitting this model to observed data:
##   y.dcasei, o.dcasei, y.dmorti, o.dmorti
calclogp <- function(params) {
    betay0 <- params[1]
    betao0 <- params[2]
    betayo0 <- params[3]
    betay1 <- params[4]
    betao1 <- params[5]
    betayo1 <- params[6]
    betay2 <- params[7]
    betao2 <- params[8]
    betayo2 <- params[9]
    ycase_rate <- params[10]
    ycase_latency <- params[11]
    yhosp_rate <- params[12]
    yhosp_latency <- params[13]
    ydied_latency <- params[14]
    ocase_rate <- params[15]
    ocase_latency <- params[16]
    ohosp_rate <- params[17]
    ohosp_latency <- params[18]
    odied_latency <- params[19]
    t0_morts <- params[20]
    t0o <- params[21]
    CLsd <- params[22]                  
    HLsd <- params[23]
    DLsd <- params[24]
    t3o <- params[25] ## start of increase, 6 June
    betay3 <- betay2
    betao3 <- betao2
    betayo3 <- betayo2
    t4o <- params[26] ## end of increase, 6 Juli
    betay4 <- params[27]
    betao4 <- params[28]
    betayo4 <- params[29]
    betay5 <- betay4 ## new lockdown, 27 Juli
    betao5 <- betao4
    betayo5 <- betayo4
    t6o <- params[30]
    betay6 <- params[31] ## end of lockdown transition 27 Juli + t6o
    betao6 <- params[32]
    betayo6 <- params[33]
    betay7 <- params[34] * betay6 ## start of school Sept 1
    betayo7 <- params[35] * betayo6
    betao7 <- params[36] * betao6

    logPriorP <- 0
    
    logPriorP <- logPriorP + dnorm(t0o, mean=0, sd=10, log=T)
    logPriorP <- logPriorP + dnorm(lockdown_offset + lockdown_transition_period + t3o, mean=d3, sd=10, log=T)
    logPriorP <- logPriorP + dnorm(lockdown_offset + lockdown_transition_period + t3o + t4o, mean=d4, sd=10, log=T)
    logPriorP <- logPriorP + dnorm(t6o, mean=G, sd=3, log=T)
    logPriorP <- logPriorP + dnorm(betay0 - betay1, mean=0, sd=2*gamma, log=T)
    logPriorP <- logPriorP + dnorm(betao0 - betao1, mean=0, sd=2*gamma, log=T)
    logPriorP <- logPriorP + dnorm(betayo0 - betayo1, mean=0, sd=2*gamma, log=T)
    logPriorP <- logPriorP + dnorm(betay1 - betay2, mean=0, sd=2*gamma, log=T)
    logPriorP <- logPriorP + dnorm(betao1 - betao2, mean=0, sd=2*gamma, log=T)
    logPriorP <- logPriorP + dnorm(betayo1 - betayo2, mean=0, sd=2*gamma, log=T)
    logPriorP <- logPriorP + dnorm(betay3 - betay4, mean=0, sd=0.5*gamma, log=T)
    logPriorP <- logPriorP + dnorm(betao3 - betao4, mean=0, sd=0.5*gamma, log=T)
    logPriorP <- logPriorP + dnorm(betayo3 - betayo4, mean=0, sd=0.5*gamma, log=T)
    logPriorP <- logPriorP + dnorm(CLsd, mean=5, sd=1, log=T)
    logPriorP <- logPriorP + dnorm(HLsd, mean=5, sd=1, log=T)
    logPriorP <- logPriorP + dnorm(DLsd, mean=5, sd=1, log=T)
    logPriorP <- logPriorP + dnorm(ydied_latency, mean=21, sd=4, log=T)
    logPriorP <- logPriorP + dnorm(odied_latency, mean=21, sd=4, log=T)

    logPriorP <- logPriorP + dnorm(betay6, mean=0.6, sd=0.2, log=T)
    logPriorP <- logPriorP + dnorm(betao6, mean=0.22, sd=0.2, log=T)
    logPriorP <- logPriorP + dnorm(betayo6, mean=0.26, sd=0.2, log=T)

    logPriorP <- logPriorP + dnorm(betay7, mean=0.85, sd=0.2, log=T)
    logPriorP <- logPriorP + dnorm(betao7, mean=0.19, sd=0.05, log=T)
    logPriorP <- logPriorP + dnorm(betayo7, mean=0.02, sd=0.02, log=T)

    logPriorP
}

calclogl <- function(params, x) {
    state <<- calculateModel(params, FitTotalPeriod)

    if (state$offset == InvalidDataOffset)
        state$offset = 1

    ## cases
    dstart <- state$offset
    dend <- state$offset + length(y.dcasei) - 1

    if (dend > length(state$y.casei)) {
        ##print("=========================== Increase FitTotalPeriod ===================")
        dend <- length(state$y.casei)
        dstart <- dend - length(y.dcasei) + 1
    }

    di <- y.dcasei

    if (dstart < 1) {
        ##print("=========================== Increase Padding ? ===================")
        dstart <- 1
        dend <- dstart + length(y.dcasei) - 1
    }

    y.loglC <- sum(dnbinom(y.dcasei[1:(d.reliable.cases-1)],
                           mu=pmax(0.1, state$y.casei[dstart:(dstart + d.reliable.cases)]),
                           size=case_nbinom_size1, log=T)) +
               sum(dnbinom(y.dcasei[d.reliable.cases:length(y.dcasei)],
                           mu=pmax(0.1, state$y.casei[(dstart + d.reliable.cases + 1):dend]),
                           size=case_nbinom_size2, log=T))

    o.loglC <- sum(dnbinom(o.dcasei[1:(d.reliable.cases-1)],
                           mu=pmax(0.1, state$o.casei[dstart:(dstart + d.reliable.cases)]),
                           size=case_nbinom_size1, log=T)) +
               sum(dnbinom(o.dcasei[d.reliable.cases:length(o.dcasei)],
                           mu=pmax(0.1, state$o.casei[(dstart + d.reliable.cases + 1):dend]),
                           size=case_nbinom_size2, log=T))

    ## hosp
    dstart <- state$offset
    dend <- state$offset + length(dhospi) - 1

    if (dend > length(state$y.hospi)) {
        ##print("=========================== Increase FitTotalPeriod ===================")
        dend <- length(state$y.hospi)
        dstart <- dend - length(dhospi) + 1
    }

    if (dstart < 1) {
        ##print("=========================== Increase Padding ? ===================")
        dstart <- 1
        dend <- dstart + length(dhospi) - 1
    }

    loglH <- sum(dnbinom(dhospi,
                         mu=pmax(0.001, (state$y.hospi + state$o.hospi)[dstart:dend]),
                         size=hosp_nbinom_size, log=T))

    ## deaths
    dstart <- state$offset
    dend <- state$offset + length(y.dmorti) - 1

    if (dend > length(state$y.deadi)) {
        ##print("=========================== Increase FitTotalPeriod ===================")
        dend <- length(state$y.deadi)
        dstart <- dend - length(y.dmorti) + 1
    }

    if (dstart < 1) {
        ##print("=========================== Increase Padding ? ===================")
        dstart <- 1
        dend <- dstart + length(y.dmorti) - 1
    }

    y.loglD <- sum(dnbinom(y.dmorti,
                           mu=pmax(0.001, state$y.deadi[dstart:dend]),
                           size=mort_nbinom_size, log=T))

    o.loglD <- sum(dnbinom(o.dmorti,
                           mu=pmax(0.001, state$o.deadi[dstart:dend]),
                           size=mort_nbinom_size, log=T))

    it <<- it + 1

    result <- y.loglC + o.loglC + loglH + y.loglD + o.loglD
    
    if (it %% 1000 == 0) {
        print(params)
	print(c(it, result))
        state <<- calcNominalState(state)
	graphs()
    }

    result
}

fit.paramnames <- c("betay0", "betao0", "betayo0",
                    "betay1", "betao1", "betayo1",
                    "betay2", "betao2", "betayo2",
                    "y.CR", "y.CL", "y.HR", "y.HL", "y.DL",
                    "o.CR", "o.CL", "o.HR", "o.HL", "o.DL",
                    "t0_morts", "t0o",
                    "CRsd", "HLsd", "DLsd",
                    "t3o",
                    "t4o", "betay4", "betao4", "betayo4",
                    "t6o", "betay6", "betao6", "betayo6",
                    "fy7", "fo7", "fyo7")
keyparamnames <- c("betay6", "betao6", "betayo6", "fy7", "fo7", "fyo7")
fitkeyparamnames <- keyparamnames

init <- c(3.6 * gamma, 3.6 * gamma, 3.6 * gamma,
          2.0 * gamma, 2.0 * gamma, 2.0 * gamma,
          0.8 * gamma, 0.8 * gamma, 0.8 * gamma,
          100, 10, 1, 10, 21,
          20, 10, 1, 10, 21,
          total_deaths_at_lockdown, -1, 5, 5, 5,
          d3 - lockdown_offset - lockdown_transition_period,
          d4 - d3, 0.8 * gamma, 0.8 * gamma, 0.8 * gamma,
          G, 0.8 * gamma, 0.8 * gamma, 0.8 * gamma,
          1.2, 1.01, 1.01)

print(init)

df_params <- data.frame(name = fit.paramnames,
                        min = c(2 * gamma, 2 * gamma, 2 * gamma,
                                1 * gamma, 1 * gamma, 1 * gamma,
                                0.05 * gamma, 0.01 * gamma, 0.01 * gamma,
                                3, 3, 0.5, 3, 10,
                                3, 3, 0.5, 3, 10,
                                0, -30, 2, 2, 2,
                                60,
                                10, 0.05 * gamma, 0.01 * gamma, 0.01 * gamma,
                                3, 0.05 * gamma, 0.01 * gamma, 0.01 * gamma,
                                1, 1, 1),
                        max = c(8 * gamma, 8 * gamma, 8 * gamma,
                                5 * gamma, 5 * gamma, 5 * gamma,
                                2 * gamma, 2 * gamma, 2 * gamma,
                                5000, 25, 2, 25, 50,
                                100, 25, 2, 25, 50,
                                max(dmort[length(dmort)] / 10, total_deaths_at_lockdown * 10),
                                30, 9, 9, 9,
                                90,
                                50, 3 * gamma, 3 * gamma, 3 * gamma,
                                12, 1.5 * gamma, 1.5 * gamma, 1.5 * gamma,
                                2, 6, 6),
                        init = init)
