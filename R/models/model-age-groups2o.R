library("deSolve")

cppname <- "model-ageo"

system(paste("R CMD SHLIB ", cppname, ".cpp", sep=""))
dyn.load(paste(cppname, ".so", sep=""))

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
    yhosp_rate <- params[11]
    yhosp_latency <- params[12]
    ydied_latency <- params[13]
    ocase_rate <- params[14]
    ohosp_rate <- params[15]
    ohosp_latency <- params[16]
    odied_latency <- params[17]
    t0_morts <- params[18]
    t0o <- params[19]
    HLsd <- params[20]
    DLsd <- params[21]
    fyifr <- params[22]
    fyhr <- params[23]
    t3o <- params[24] ## start of increase, ~ 6 June
    betay3 <- betay2
    t4o <- params[25] ## end of increase, ~ 6 July
    betay4 <- params[26]
    betao4 <- params[27]
    betayo4 <- params[28]
    betao3 <- betao4
    betayo3 <- betayo4
    betay5 <- params[29] ## new lockdown, 27 July
    betao5 <- params[30]
    betayo5 <- params[31]
    t6o <- params[32] ## exiting lockdown ~ mid August
    betay6 <- params[33]
    betao6 <- params[34]
    betayo6 <- params[35]
    t7o <- params[36]  ## sept 15 ?
    betay7 <- params[37]
    betao7 <- params[38]
    betayo7 <- params[39]
    betay8 <- betay7
    betao8 <- betao7
    betayo8 <- betayo7
    betay9 <- params[40] * betay8
    betao9 <- params[41] * betao8
    betayo9 <- params[42] * betayo8
    ifrred <- params[43]
    ycase_latency <- params[44]
    ocase_latency <- params[45]

    betay10 <- betay9
    betao10 <- betao9
    betayo10 <- betayo9

    ## 2nd lockdown : 1.2 ~ 1st lockdown
    betay11 <- 1.2 * betay2
    betao11 <- min(1.2 * betao2, betao9)
    betayo11 <- min(1.2 * betayo2, betayo9)
    
    CLsd <- 5
    
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
               t7 = 1E10, betay7 = betay7, betao7 = betao7, betayo7 = betayo7,
               t8 = 1E10, betay8 = betay8, betao8 = betao8, betayo8 = betayo8,
               t9 = 1E10, betay9 = betay9, betao9 = betao9, betayo9 = betayo9,
               t10 = 1E10, betay10 = betay10, betao10 = betao10, betayo10 = betayo10,
               t11 = 1E10, betay11 = betay11, betao11 = betao11, betayo11 = betayo11)

    Y <- c(Sy = y.N - Initial, E1y = Initial, E2y = 0, Iy = 0, Ry = 0,
           So = o.N, E1o = 0, E2o = 0, Io = 0, Ro = 0)

    period1 <- 60
    times <- (padding + 1):(padding + period1)

    out <- ode(Y, times, func = "derivs", parms = parms,
               dllname = cppname,
               initfunc = "initmod", nout = 4, outnames = c("Re", "Rt", "y.Re", "o.Re"))

    state$y.S[(padding + 1):(padding + period1)] = out[,2]
    state$o.S[(padding + 1):(padding + period1)] = out[,7]

    s2 <- convolute(state$y.S, padding + 1, padding + period1, y.died_cv_profile)
    state$y.died[(padding + 1):(padding + period1)] = (y.N - s2) * y.died_rate

    s2 <- convolute(state$o.S, padding + 1, padding + period1, o.died_cv_profile)
    state$o.died[(padding + 1):(padding + period1)] = (o.N - s2) * o.died_rate

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
        t7 <- t6 + t7o
        t8 <- data_offset + d8
        t9 <- t8 + 7
        t10 <- data_offset + d10
        t11 <- data_offset + d10 + 4

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
                   t7 = t7, betay7 = betay7, betao7 = betao7, betayo7 = betayo7,
                   t8 = t8, betay8 = betay8, betao8 = betao8, betayo8 = betayo8,
                   t9 = t9, betay9 = betay9, betao9 = betao9, betayo9 = betayo9,
                   t10 = t10, betay10 = betay10, betao10 = betao10, betayo10 = betayo10,
                   t11 = t11, betay11 = betay11, betao11 = betao11, betayo11 = betayo11)

        times <- (padding + 1):(padding + period)
        
        out <- ode(Y, times, func = "derivs", parms = parms,
                   dllname = cppname,
                   initfunc = "initmod", nout = 4, outnames = c("Re", "Rt", "y.Re", "o.Re"))
    }

    state$y.S[(padding + 1):(padding + period)] = out[,2]
    state$y.E[(padding + 1):(padding + period)] = out[,3] + out[,4]
    state$y.I[(padding + 1):(padding + period)] = out[,5]
    state$y.R[(padding + 1):(padding + period)] = out[,6]
    state$o.S[(padding + 1):(padding + period)] = out[,7]
    state$o.E[(padding + 1):(padding + period)] = out[,8] + out[,9]
    state$o.I[(padding + 1):(padding + period)] = out[,10]
    state$o.R[(padding + 1):(padding + period)] = out[,11]
    state$Re[(padding + 1):(padding + period)] = out[,12]
    state$Rt[(padding + 1):(padding + period)] = out[,13]
    state$y.Re[(padding + 1):(padding + period)] = out[,14]
    state$o.Re[(padding + 1):(padding + period)] = out[,15]

    gycr <- pmin(1, ycase_rate * y.ifr * (1 + (fyifr - 1) * g))     ## fyifr: age reduction
    gyifr <- y.ifr * (1 + (fyifr - 1) * g) * (1 - ifrred * gtrimp)  ## ifrred: treatment improvement -> on lowers IFR 
    gyhr <- yhosp_rate * y.hr * (1 + (fyhr - 1) * g)

    gocr <- pmin(1, ocase_rate * o.ifr)
    goifr <- o.ifr                         * (1 - ifrred * gtrimp)
    gohr <- ohosp_rate * o.hr

    t.symp <- data_offset + d.symp.cases

    ## y.case
    s2 <- y.N - convolute(state$y.S, padding + 1, padding + period, y.case_cv_profile)
    s2i <- c(s2[1], diff(s2))
    if (data_offset != InvalidDataOffset) {
        ## IFR is based on death, IHR is shifted in time (DL - CL) earlier
        t1 <- data_offset + round(-ydied_latency + ycase_latency)
        t2 <- min(padding + period, t1 + length(gyifr) - 1)
        if (t1 > padding && t2 > t1) {
            state$y.casei[(padding + 1):t1] = s2i[1:(t1 - padding)] * gycr[1]
            state$y.casei[t1:t2] = s2i[(t1 - padding):(t2 - padding)] * gycr[1:(t2 - t1 + 1)]
            state$y.casei[t2:(padding + period)] = s2i[(t2 - padding):period] * gycr[length(gycr)]
        }

        state$y.casei[t.symp:length(state$y.casei)] =
            state$y.casei[t.symp:length(state$y.casei)] * symp.cases.factor
    } else {
        state$y.casei[(padding + 1):(padding + period)] = s2i * gycr[1]
    }

    state$y.case = cumsum(state$y.casei)

    ## y.hosp
    s2 <- y.N - convolute(state$y.S, padding + 1, padding + period, y.hosp_cv_profile)
    s2i <- c(s2[1], diff(s2))
    if (data_offset != InvalidDataOffset) {
        ## IFR is based on death, IHR is shifted in time (DL - HL) earlier
        t1 <- data_offset + round(-ydied_latency + yhosp_latency)
        t2 <- min(padding + period, t1 + length(gyhr) - 1)
        if (t1 > padding && t2 > t1) {
            state$y.hospi[(padding + 1):t1] = s2i[1:(t1 - padding)] * gyhr[1]
            state$y.hospi[t1:t2] = s2i[(t1 - padding):(t2 - padding)] * gyhr[1:(t2 - t1 + 1)]
            state$y.hospi[t2:(padding + period)] = s2i[(t2 - padding):period] * gyhr[length(gyhr)]
        }
    } else {
        state$y.hospi[(padding + 1):(padding + period)] = s2i * gyhr[1]
    }

    ## y.dead
    s2 <- y.N - convolute(state$y.S, padding + 1, padding + period, y.died_cv_profile)
    s2i <- c(s2[1], diff(s2))
    if (data_offset != InvalidDataOffset) {
        t1 <- data_offset
        t2 <- min(padding + period, t1 + length(gyifr) - 1)
        if (t1 > padding && t2 > t1) {
            state$y.deadi[(padding + 1):t1] = s2i[1:(t1 - padding)] * gyifr[1]
            state$y.deadi[t1:t2] = s2i[(t1 - padding):(t2 - padding)] * gyifr[1:(t2 - t1 + 1)]
            state$y.deadi[t2:(padding + period)] = s2i[(t2 - padding):period] * gyifr[length(gyifr)]
        }
    } else {
        state$y.deadi[(padding + 1):(padding + period)] = s2i * gyifr[1]
    }
    state$y.died = cumsum(state$y.deadi)

    ## o.case
    s2 <- o.N - convolute(state$o.S, padding + 1, padding + period, o.case_cv_profile)
    s2i <- c(s2[1], diff(s2))
    if (data_offset != InvalidDataOffset) {
        t1 <- data_offset + round(-ydied_latency + ocase_latency)
        t2 <- min(padding + period, t1 + length(gocr) - 1)
        if (t1 > padding && t2 > t1) {
            state$o.casei[(padding + 1):t1] = s2i[1:(t1 - padding)] * gocr[1]
            state$o.casei[t1:t2] = s2i[(t1 - padding):(t2 - padding)] * gocr[1:(t2 - t1 + 1)]
            state$o.casei[t2:(padding + period)] = s2i[(t2 - padding):period] * gocr[length(gocr)]
        }

        state$o.casei[t.symp:length(state$o.casei)] =
            state$o.casei[t.symp:length(state$o.casei)] * symp.cases.factor

    } else {
        state$o.casei[(padding + 1):(padding + period)] = s2i * gocr[1]
    }

    state$o.case = cumsum(state$o.casei)

    ## o.hosp
    s2 <- o.N - convolute(state$o.S, padding + 1, padding + period, o.hosp_cv_profile)
    s2i <- c(s2[1], diff(s2))
    if (data_offset != InvalidDataOffset) {
        t1 <- data_offset + round(-ydied_latency + ohosp_latency)
        t2 <- min(padding + period, t1 + length(gohr) - 1)
        if (t1 > padding && t2 > t1) {
            state$o.hospi[(padding + 1):t1] = s2i[1:(t1 - padding)] * gohr[1]
            state$o.hospi[t1:t2] = s2i[(t1 - padding):(t2 - padding)] * gohr[1:(t2 - t1 + 1)]
            state$o.hospi[t2:(padding + period)] = s2i[(t2 - padding):period] * gohr[length(gohr)]
        }
    } else {
        state$o.hospi[(padding + 1):(padding + period)] = s2i * gohr[1]
    }

    ## o.dead
    s2 <- o.N - convolute(state$o.S, padding + 1, padding + period, o.died_cv_profile)
    s2i <- c(s2[1], diff(s2))
    if (data_offset != InvalidDataOffset) {
        t1 <- data_offset
        t2 <- min(padding + period, t1 + length(goifr) - 1)
        if (t1 > padding && t2 > t1) {
            state$o.deadi[(padding + 1):t1] = s2i[1:(t1 - padding)] * goifr[1]
            state$o.deadi[t1:t2] = s2i[(t1 - padding):(t2 - padding)] * goifr[1:(t2 - t1 + 1)]
            state$o.deadi[t2:(padding + period)] = s2i[(t2 - padding):period] * goifr[length(goifr)]
        }
    } else {
        state$o.deadi[(padding + 1):(padding + period)] = s2i * goifr[1]
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
    yhosp_rate <- params[11]
    yhosp_latency <- params[12]
    ydied_latency <- params[13]
    ocase_rate <- params[14]
    ohosp_rate <- params[15]
    ohosp_latency <- params[16]
    odied_latency <- params[17]
    t0_morts <- params[18]
    t0o <- params[19]
    HLsd <- params[20]
    DLsd <- params[21]
    fyifr <- params[22]
    fyhr <- params[23]
    t3o <- params[24] ## start of increase, ~ 6 June
    betay3 <- betay2
    t4o <- params[25] ## end of increase, ~ 6 July
    betay4 <- params[26]
    betao4 <- params[27]
    betayo4 <- params[28]
    betao3 <- betao4
    betayo3 <- betayo4
    betay5 <- params[29] ## new lockdown, 27 July
    betao5 <- params[30]
    betayo5 <- params[31]
    t6o <- params[32] ## exiting lockdown ~ mid August
    betay6 <- params[33]
    betao6 <- params[34]
    betayo6 <- params[35]
    t7o <- params[36]  ## sept 15 ?
    betay7 <- params[37]
    betao7 <- params[38]
    betayo7 <- params[39]
    betay8 <- betay7
    betao8 <- betao7
    betayo8 <- betayo7
    fy9 <- params[40]
    fo9 <- params[41]
    fyo9 <- params[42]
    ifrred <- params[43]
    ycase_latency <- params[44]
    ocase_latency <- params[45]

    logPriorP <- 0
    
    logPriorP <- logPriorP + dnorm(t0o, mean=0, sd=10, log=T)
    logPriorP <- logPriorP + dnorm(lockdown_offset + lockdown_transition_period + t3o,
                                   mean=d3, sd=10, log=T)
    logPriorP <- logPriorP + dnorm(lockdown_offset + lockdown_transition_period + t3o + t4o,
                                   mean=d4, sd=10, log=T)
    logPriorP <- logPriorP + dnorm(d5 + t6o, mean=d6, sd=5, log=T)
    logPriorP <- logPriorP + dnorm(d5 + t6o + t7o, mean=d7, sd=4, log=T)

    SD <- log(2) ## SD = */ 2
    lSD <- log(1.5) ## SD = */ 1.5

    logPriorP <- logPriorP + dlnorm(betay0/betay1, 0, SD, log=T)
    logPriorP <- logPriorP + dlnorm(betao0/betao1, 0, SD, log=T)
    logPriorP <- logPriorP + dlnorm(betayo0/betayo1, 0, SD, log=T)
    logPriorP <- logPriorP + dlnorm(betay1/betay2, 0, SD, log=T)
    logPriorP <- logPriorP + dlnorm(betao1/betao2, 0, SD, log=T)
    logPriorP <- logPriorP + dlnorm(betayo1/betayo2, 0, SD, log=T)
    logPriorP <- logPriorP + dlnorm(betay3/betay4, 0, SD, log=T)
    logPriorP <- logPriorP + dlnorm(betao2/betao4, 0, SD, log=T)
    logPriorP <- logPriorP + dlnorm(betayo2/betayo4, 0, SD, log=T)
    logPriorP <- logPriorP + dlnorm(betay5/betay6, 0, SD, log=T)
    logPriorP <- logPriorP + dlnorm(betao5/betao6, 0, SD, log=T)
    logPriorP <- logPriorP + dlnorm(betayo5/betayo6, 0, SD, log=T)
    logPriorP <- logPriorP + dlnorm(betay6/betay7, 0, lSD, log=T)
    logPriorP <- logPriorP + dlnorm(betao6/betao7, 0, lSD, log=T)
    logPriorP <- logPriorP + dlnorm(betayo6/betayo7, 0, lSD, log=T)

    logPriorP <- logPriorP + dnorm(betao7, mean=0.5 * gamma, sd=0.2 * gamma, log=T)
    logPriorP <- logPriorP + dnorm(betao8 * fo9, mean=0.5 * gamma, sd=0.2 * gamma, log=T)
    
    logPriorP <- logPriorP + dlnorm(fy9, log(0.85), log(1.07), log=T)
    logPriorP <- logPriorP + dlnorm(fo9, log(0.85), log(1.07), log=T)
    logPriorP <- logPriorP + dlnorm(fyo9, log(0.85), log(1.07), log=T)

    logPriorP <- logPriorP + dnorm(HLsd, mean=5, sd=1, log=T)
    logPriorP <- logPriorP + dnorm(DLsd, mean=5, sd=1, log=T)
    logPriorP <- logPriorP + dnorm(ycase_latency, mean=10, sd=2, log=T)
    logPriorP <- logPriorP + dnorm(ocase_latency, mean=10, sd=2, log=T)
    logPriorP <- logPriorP + dnorm(yhosp_latency, mean=13, sd=3, log=T)
    logPriorP <- logPriorP + dnorm(ohosp_latency, mean=13, sd=3, log=T)
    logPriorP <- logPriorP + dnorm(ydied_latency, mean=21, sd=0.5, log=T)
    logPriorP <- logPriorP + dnorm(odied_latency, mean=21, sd=0.5, log=T)

    logPriorP <- logPriorP + dlnorm(fyifr, log(0.3), log(2), log=T) # */ 2
    logPriorP <- logPriorP + dlnorm(fyhr, log(0.3), log(2), log=T)
    logPriorP <- logPriorP + dlnorm(ifrred, log(0.35), log(1.3), log=T) # */ 1.3

    logPriorP <- logPriorP + dlnorm(yhosp_rate, 0, lSD, log=T)
    logPriorP <- logPriorP + dlnorm(ohosp_rate, 0, lSD, log=T)
    
    logPriorP
}

calclogl <- function(params, x) {
    it <<- it + 1

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
                           size=case_nbinom_sizey2, log=T))

    o.loglC <- sum(dnbinom(o.dcasei[1:(d.reliable.cases-1)],
                           mu=pmax(0.1, state$o.casei[dstart:(dstart + d.reliable.cases)]),
                           size=case_nbinom_size1, log=T)) +
               sum(dnbinom(o.dcasei[d.reliable.cases:length(o.dcasei)],
                           mu=pmax(0.1, state$o.casei[(dstart + d.reliable.cases + 1):dend]),
                           size=case_nbinom_sizeo2, log=T))

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

    loglH <- sum(dnbinom(dhospi[1:(d.reliable.hosp-1)],
                         mu=pmax(0.1, (state$y.hospi + state$o.hospi)
                                 [dstart:(dstart + d.reliable.hosp)]),
                         size=hosp_nbinom_size1, log=T)) +
             sum(dnbinom(dhospi[d.reliable.hosp:length(dhosp)],
                         mu=pmax(0.1, (state$y.hospi + state$o.hospi)
                                 [(dstart + d.reliable.hosp + 1):dend]),
                         size=hosp_nbinom_size2, log=T)) +
             sum(dnbinom(dhospi[(length(dhosp) - 7):length(dhosp)],
                         mu=pmax(0.1, (state$y.hospi + state$o.hospi)[(dend - 7):dend]),
                         size=hosp_nbinom_size2 * hosp_last7, log=T))

    if (it %% 1000 == 0) {
        a <- dhospi[length(dhosp)]
        b <- (state$y.hospi + state$o.hospi)[dend]
        print(c(a, b, dnbinom(a, mu=b, size=hosp_nbinom_size2 * hosp_last7, log=T)))
    }
    
    ## between June 22st and September 14th, ratio y/o should be 56.7/43.3 +/- 5%
    ytotal = sum(state$y.hospi[(dstart + d.hosp.o1):(dstart + d.hosp.o2)])
    ototal = sum(state$o.hospi[(dstart + d.hosp.o1):(dstart + d.hosp.o2)])
    loglHRatio <- dlnorm(ytotal/(ytotal + ototal), log(56.7/100), log(1.05), log=T)
    ##print(c(loglH, ytotal, ototal, loglHRatio))
    
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
                           size=mort_nbinom_size*2, log=T))

    o.loglD <- sum(dnbinom(o.dmorti,
                           mu=pmax(0.001, state$o.deadi[dstart:dend]),
                           size=o.wdmorti * mort_nbinom_size, log=T))

    result <- y.loglC + o.loglC + loglH + loglHRatio + y.loglD + o.loglD
    
    if (it %% 1000 == 0) {
        print(params)
	print(c(it, y.loglC, o.loglC, loglH, loglHRatio, y.loglD, o.loglD, result))
        state <<- calcNominalState(state)
	graphs()
    }

    result
}

fit.paramnames <- c("betay0", "betao0", "betayo0",
                    "betay1", "betao1", "betayo1",
                    "betay2", "betao2", "betayo2",
                    "y.CR", "y.HR", "y.HL", "y.DL",
                    "o.CR", "o.HR", "o.HL", "o.DL",
                    "t0_morts", "t0o",
                    "HLsd", "DLsd",
                    "fyifr", "fyhr",
                    "t3o",
                    "t4o", "betay4", "betao4", "betayo4",
                    "betay5", "betao5", "betayo5",
                    "t6o", "betay6", "betao6", "betayo6",
                    "t7o", "betay7", "betao7", "betayo7",
                    "fy9", "fo9", "fyo9",
                    "ifrred", "y.CL", "o.CL")
keyparamnames <- c("betay6", "betao6", "betayo6",
                   "betay7", "betao7", "betayo7",
                   "fy9", "fo9", "fyo9", "ifrred")
fitkeyparamnames <- keyparamnames

init <- c(2.9, 1.5, 1.5,
          1.5, 0.8, 0.7,
          0.7, 0.5, 0.07,
          2200, 1.2, 15, 19,
          60, 1.3, 13, 21,
          14, -8, 4.1, 6, 0.13, 0.10,
          d3 - lockdown_offset - lockdown_transition_period,
          d4 - d3, 1, 0.11, 0.04,
                   0.6, 0.16, 0.03,
          d6 - d5, 1.2, 0.3, 0.08,
          d7 - d6, 1.7, 0.4, 0.12,
          0.85, 0.85, 0.85, 0.3, 11, 10)

df_params <- data.frame(name = fit.paramnames,
                        min = c(2 * gamma, 1 * gamma, 1 * gamma,
                                1 * gamma, 0.5 * gamma, 0.5 * gamma,
                                0.1 * gamma, 0.1 * gamma, 0.01 * gamma,
                                500, 0.5, 6, 14,
                                3, 0.5, 6, 14,
                                0, -20, 2, 2, 0.05, 0.05,
                                60,
                                10, 0.5 * gamma, 0.01 * gamma, 0.002 * gamma,
                                    0.1 * gamma, 0.01 * gamma, 0.002 * gamma,
                                10, 0.5 * gamma, 0.01 * gamma, 0.002 * gamma,
                                10, 0.5 * gamma, 0.01 * gamma, 0.002 * gamma,
                                0.5, 0.5, 0.5, 0.1, 4, 4),
                        max = c(8 * gamma, 5 * gamma, 5 * gamma,
                                5 * gamma, 3 * gamma, 3 * gamma,
                                2 * gamma, 2 * gamma, 2 * gamma,
                                5000, 2, 20, 25,
                                100, 6, 20, 25,
                                max(dmort[length(dmort)] / 10, total_deaths_at_lockdown * 10),
                                10, 9, 9, 1, 1,
                                100,
                                50, 3 * gamma, 1 * gamma, 0.5 * gamma,
                                    2 * gamma, 1 * gamma, 0.5 * gamma,
                                50, 3 * gamma, 3 * gamma, 0.5 * gamma,
                                50, 3 * gamma, 3 * gamma, 0.5 * gamma,
                                1.2, 1.2, 1.2, 0.7, 14, 17),
                        init = init)

print(df_params)

