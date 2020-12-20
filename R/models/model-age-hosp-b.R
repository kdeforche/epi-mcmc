library("deSolve")

cppname <- "model-ager"

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
eta <- 1/(0.75 * 356) ## 9 months
CLsd <- 2.5

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
    betay3 <- params[25]
    betao3 <- params[26]
    betayo3 <- params[27]
    t4o <- params[28] ## end of increase, ~ 6 July
    betay4 <- params[29]
    betao4 <- params[30]
    betayo4 <- params[31]
    betay5 <- params[32] ## new lockdown, 27 July
    betao5 <- params[33]
    betayo5 <- params[34]
    t6o <- params[35] ## exiting lockdown ~ mid August
    betay6 <- params[36]
    betao6 <- params[37]
    betayo6 <- params[38]
    t7o <- params[39]  ## sept 15 ?
    betay7 <- params[40]
    betao7 <- params[41]
    betayo7 <- params[42]
    betay8 <- betay7
    betao8 <- betao7
    betayo8 <- betayo7
    betay9 <- params[43]
    betao9 <- params[44]
    betayo9 <- params[45]
    ifrred <- params[46]
    ycase_latency <- params[47]
    ocase_latency <- params[48]
    t10o <- params[49]
    t11o <- params[50]
    betay11 <- params[51]
    betao11 <- params[52]
    betayo11 <- params[53]

    betay10 <- betay9
    betao10 <- betao9
    betayo10 <- betayo9

    betay12 <- betay9
    betao12 <- betao9
    betayo12 <- betayo9
    
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
    state$y.i <- rep(0, padding + period)
    state$o.i <- rep(0, padding + period)
    state$y.beta <- rep(0, padding + period)
    state$o.beta <- rep(0, padding + period)
    state$yo.beta <- rep(0, padding + period)

    parms <- c(Ny = y.N, No = o.N,
               a1 = a1, a2 = a2, gamma = gamma, eta = eta,
               tuncertain = 1E10, funcertain=1,
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
               t11 = 1E10, betay11 = betay11, betao11 = betao11, betayo11 = betayo11,
               t12 = 1E10, betay12 = betay12, betao12 = betao12, betayo12 = betayo12)

    Y <- c(Sy = y.N - Initial, E1y = Initial, E2y = 0, Iy = 0, Ry = 0,
           So = o.N, E1o = 0, E2o = 0, Io = 0, Ro = 0)

    period1 <- 60
    times <- (padding + 1):(padding + period1)

    out <- ode(Y, times, func = "derivs", parms = parms,
               dllname = cppname,
               initfunc = "initmod", nout = 9,
               outnames = c("Re", "Rt", "y.Re", "o.Re", "y.i", "o.i",
                            "y.beta", "o.beta", "yo.beta"))

    state$y.i[(padding + 1):(padding + period1)] = out[,16]
    state$o.i[(padding + 1):(padding + period1)] = out[,17]

    y.deadi <- convolute(state$y.i, padding + 1, padding + period1, y.died_cv_profile) * y.died_rate
    state$y.died[(padding + 1):(padding + period1)] = cumsum(y.deadi)

    o.deadi <- convolute(state$o.i, padding + 1, padding + period1, o.died_cv_profile) * o.died_rate
    state$o.died[(padding + 1):(padding + period1)] = cumsum(o.deadi)

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
        t9 <- data_offset + d9
        t10 <- data_offset + t10o
        t11 <- data_offset + t11o
        t12 <- data_offset + d12
        tuncertain <- data_offset + duncertain
        funcertain <- 1 ## rlnorm(1, meanlog=0, sdlog=log(1.1))

        parms <- c(Ny = y.N, No = o.N,
                   a1 = a1, a2 = a2, gamma = gamma, eta = eta,
                   tuncertain=tuncertain, funcertain=funcertain,
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
                   t11 = t11, betay11 = betay11, betao11 = betao11, betayo11 = betayo11,
                   t12 = t12, betay12 = betay12, betao12 = betao12, betayo12 = betayo12)

        times <- (padding + 1):(padding + period)
        
        out <- ode(Y, times, func = "derivs", parms = parms,
                   dllname = cppname,
                   initfunc = "initmod", nout = 9,
                   outnames = c("Re", "Rt", "y.Re", "o.Re", "y.i", "o.i",
                                "y.beta", "o.beta", "yo.beta"))
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
    state$y.i[(padding + 1):(padding + period)] = out[,16]
    state$o.i[(padding + 1):(padding + period)] = out[,17]
    state$y.beta[(padding + 1):(padding + period)] = out[,18]
    state$o.beta[(padding + 1):(padding + period)] = out[,19]
    state$yo.beta[(padding + 1):(padding + period)] = out[,20]

    gycr <- pmin(1, ycase_rate * y.ifr * (1 + (fyifr - 1) * g))     ## fyifr: age reduction
    gyifr <- y.ifr * (1 + (fyifr - 1) * g) * (1 - ifrred * gtrimp)  ## ifrred: treatment improvement -> on lowers IFR 
    gyhr <- yhosp_rate * y.hr * (1 + (fyhr - 1) * g)

    gocr <- pmin(1, ocase_rate * o.ifr)
    goifr <- o.ifr                         * (1 - ifrred * gtrimp)
    gohr <- ohosp_rate * o.hr

    t.symp <- data_offset + d.symp.cases
    t.all <- data_offset + d.all.cases - 1
    
    ## y.case
    s2i <- convolute(state$y.i, padding + 1, padding + period, y.case_cv_profile)
    if (data_offset != InvalidDataOffset) {
        ## IFR is based on death, IHR is shifted in time (DL - CL) earlier
        t1 <- data_offset + round(-ydied_latency + ycase_latency)
        t2 <- min(padding + period, t1 + length(gyifr) - 1)
        if (t1 > padding && t2 > t1) {
            state$y.casei[(padding + 1):t1] = s2i[1:(t1 - padding)] * gycr[1]
            state$y.casei[t1:t2] = s2i[(t1 - padding):(t2 - padding)] * gycr[1:(t2 - t1 + 1)]
            state$y.casei[t2:(padding + period)] = s2i[(t2 - padding):period] * gycr[length(gycr)]
        }

        state$y.casei[t.symp:min(t.all, length(state$y.casei))] =
            state$y.casei[t.symp:min(t.all, length(state$y.casei))] * y.symp.profile
    } else {
        state$y.casei[(padding + 1):(padding + period)] = s2i * gycr[1]
    }
    state$y.case = cumsum(state$y.casei)

    ## y.hosp
    s2i <- convolute(state$y.i, padding + 1, padding + period, y.hosp_cv_profile)
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
    s2i <- convolute(state$y.i, padding + 1, padding + period, y.died_cv_profile)
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
    s2i <- convolute(state$o.i, padding + 1, padding + period, o.case_cv_profile)
    if (data_offset != InvalidDataOffset) {
        t1 <- data_offset + round(-ydied_latency + ocase_latency)
        t2 <- min(padding + period, t1 + length(gocr) - 1)
        if (t1 > padding && t2 > t1) {
            state$o.casei[(padding + 1):t1] = s2i[1:(t1 - padding)] * gocr[1]
            state$o.casei[t1:t2] = s2i[(t1 - padding):(t2 - padding)] * gocr[1:(t2 - t1 + 1)]
            state$o.casei[t2:(padding + period)] = s2i[(t2 - padding):period] * gocr[length(gocr)]
        }
        state$o.casei[t.symp:min(t.all, length(state$o.casei))] =
            state$o.casei[t.symp:min(t.all, length(state$o.casei))] * o.symp.profile
    } else {
        state$o.casei[(padding + 1):(padding + period)] = s2i * gocr[1]
    }
    state$o.case = cumsum(state$o.casei)

    ## o.hosp
    s2i <- convolute(state$o.i, padding + 1, padding + period, o.hosp_cv_profile)
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
    s2i <- convolute(state$o.i, padding + 1, padding + period, o.died_cv_profile)
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

    posterior$y.ld2 = posterior$betay11 / posterior$betay9
    posterior$o.ld2 = posterior$betao11 / posterior$betao9
    posterior$yo.ld2 = posterior$betayo11 / posterior$betayo9
    
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
    betay3 <- params[25]
    betao3 <- params[26]
    betayo3 <- params[27]
    t4o <- params[28] ## end of increase, ~ 6 July
    betay4 <- params[29]
    betao4 <- params[30]
    betayo4 <- params[31]
    betay5 <- params[32] ## new lockdown, 27 July
    betao5 <- params[33]
    betayo5 <- params[34]
    t6o <- params[35] ## exiting lockdown ~ mid August
    betay6 <- params[36]
    betao6 <- params[37]
    betayo6 <- params[38]
    t7o <- params[39]  ## sept 15 ?
    betay7 <- params[40]
    betao7 <- params[41]
    betayo7 <- params[42]
    betay8 <- betay7
    betao8 <- betao7
    betayo8 <- betayo7
    betay9 <- params[43]
    betao9 <- params[44]
    betayo9 <- params[45]
    ifrred <- params[46]
    ycase_latency <- params[47]
    ocase_latency <- params[48]
    t10o <- params[49]
    t11o <- params[50]
    betay11 <- params[51]
    betao11 <- params[52]
    betayo11 <- params[53]

    logPriorP <- 0

    ## easier is to have all offsets relative to startdate ?
    logPriorP <- logPriorP + dnorm(t0o, mean=0, sd=10, log=T)
    logPriorP <- logPriorP + dnorm(lockdown_offset + lockdown_transition_period + t3o,
                                   mean=d3, sd=10, log=T)
    logPriorP <- logPriorP + dnorm(lockdown_offset + lockdown_transition_period + t3o + t4o,
                                   mean=d4, sd=10, log=T)
    logPriorP <- logPriorP + dnorm(d5 + t6o, mean=d6, sd=5, log=T)
    logPriorP <- logPriorP + dnorm(d5 + t6o + t7o, mean=d7, sd=2, log=T)
    logPriorP <- logPriorP + dnorm(t10o, mean=d10, sd=7, log=T)
    logPriorP <- logPriorP + dnorm(t11o, mean=d11, sd=7, log=T)

    SD <- log(2) ## SD = */ 2
    lSD <- log(1.3) ## SD = */ 1.3
    llSD <- log(1.2) ## SD = */ 1.2

    logPriorP <- logPriorP + dlnorm(betay0/betay1, 0, llSD, log=T)
    logPriorP <- logPriorP + dlnorm(betao0/betao1, 0, llSD, log=T)
    logPriorP <- logPriorP + dlnorm(betayo0/betayo1, 0, llSD, log=T)
    logPriorP <- logPriorP + dlnorm(betay1/betay2, 0, SD, log=T)
    logPriorP <- logPriorP + dlnorm(betao1/betao2, 0, SD, log=T)
    logPriorP <- logPriorP + dlnorm(betayo1/betayo2, 0, SD, log=T)
    logPriorP <- logPriorP + dlnorm(betay2/betay3, 0, SD, log=T)
    logPriorP <- logPriorP + dlnorm(betao2/betao3, 0, SD, log=T)
    logPriorP <- logPriorP + dlnorm(betayo2/betayo3, 0, SD, log=T)
    logPriorP <- logPriorP + dlnorm(betay3/betay4, 0, SD, log=T)
    logPriorP <- logPriorP + dlnorm(betao3/betao4, 0, SD, log=T)
    logPriorP <- logPriorP + dlnorm(betayo3/betayo4, 0, SD, log=T)
    logPriorP <- logPriorP + dlnorm(betay5/betay6, 0, SD, log=T)
    logPriorP <- logPriorP + dlnorm(betao5/betao6, 0, SD, log=T)
    logPriorP <- logPriorP + dlnorm(betayo5/betayo6, 0, SD, log=T)
    logPriorP <- logPriorP + dlnorm(betay6/betay7, 0, SD, log=T)
    logPriorP <- logPriorP + dlnorm(betao6/betao7, 0, SD, log=T)
    logPriorP <- logPriorP + dlnorm(betayo6/betayo7, 0, SD, log=T)
    logPriorP <- logPriorP + dlnorm(betay7/betay9, log(1.3), lSD, log=T)
    logPriorP <- logPriorP + dlnorm(betao7/betao9, log(1.3), lSD, log=T)
    logPriorP <- logPriorP + dlnorm(betayo7/betayo9, log(1.3), lSD, log=T)

    logPriorP <- logPriorP + dlnorm(betay9/betay11, 0, SD, log=T)
    logPriorP <- logPriorP + dlnorm(betao9/betao11, 0, SD, log=T)
    logPriorP <- logPriorP + dlnorm(betayo9/betayo11, 0, SD, log=T)
    
    logPriorP <- logPriorP + dnorm(HLsd, mean=4, sd=1, log=T)
    logPriorP <- logPriorP + dnorm(DLsd, mean=5, sd=1, log=T)
    logPriorP <- logPriorP + dnorm(ycase_latency, mean=7, sd=3, log=T)
    logPriorP <- logPriorP + dnorm(ocase_latency, mean=7, sd=3, log=T)
    logPriorP <- logPriorP + dnorm(yhosp_latency, mean=13, sd=1, log=T)
    logPriorP <- logPriorP + dnorm(ohosp_latency, mean=13, sd=1, log=T)
    logPriorP <- logPriorP + dnorm(ydied_latency, mean=21, sd=0.5, log=T)
    logPriorP <- logPriorP + dnorm(odied_latency, mean=21, sd=0.5, log=T)

    ## .2
    logPriorP <- logPriorP + dlnorm(fyifr, log(0.6), log(1.5), log=T) # */ 1.2
    ##logPriorP <- logPriorP + dlnorm(fyifr, log(0.35), log(1.2), log=T) # */ 2

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

    last <- 21
    
    y.loglC <- sum(dnbinom(y.dcasei[1:(d.reliable.cases-1)],
                           mu=pmax(0.1, state$y.casei[dstart:(dstart + d.reliable.cases)]),
                           size=case_nbinom_size1, log=T)) +
               sum(dnbinom(y.dcasei[d.reliable.cases:length(y.dcasei)],
                           mu=pmax(0.1, state$y.casei[(dstart + d.reliable.cases + 1):dend]),
                           size=case_nbinom_sizey2, log=T)) +
               sum(dnbinom(y.dcasei[(length(y.dcasei) - last):length(y.dcasei)],
                           mu=pmax(0.1, state$y.casei[(dend - last):dend]),
                           size=case_nbinom_sizey2 * f_last7, log=T))

    o.loglC <- sum(dnbinom(o.dcasei[1:(d.reliable.cases-1)],
                           mu=pmax(0.1, state$o.casei[dstart:(dstart + d.reliable.cases)]),
                           size=case_nbinom_size1, log=T)) +
               sum(dnbinom(o.dcasei[d.reliable.cases:length(o.dcasei)],
                           mu=pmax(0.1, state$o.casei[(dstart + d.reliable.cases + 1):dend]),
                           size=case_nbinom_sizeo2, log=T)) +
               sum(dnbinom(o.dcasei[(length(o.dcasei) - last):length(o.dcasei)],
                           mu=pmax(0.1, state$o.casei[(dend - last):dend]),
                           size=case_nbinom_sizeo2 * f_last7, log=T))

    ## hosp
    dstart <- state$offset
    dend <- state$offset + length(ws) - 1

    if (dend > length(state$y.hospi)) {
        ##print("=========================== Increase FitTotalPeriod ===================")
        dend <- length(state$y.hospi)
        dstart <- dend - length(ws) + 1
    }

    if (dstart < 1) {
        ##print("=========================== Increase Padding ? ===================")
        dstart <- 1
        dend <- dstart + length(ws) - 1
    }

    ## state dstart-dend corresponds with dhospi

    state.y.whospi <- aggregate(state$y.hospi[dstart:dend], by=list(week=ws), FUN=sum)
    state.o.whospi <- aggregate(state$o.hospi[dstart:dend], by=list(week=ws), FUN=sum)
    
    loglH.1 <- sum(dnbinom(y.whospi, mu=pmax(0.1, state.y.whospi$x),
                           size=whosp_nbinom_size, log=T)) +
        sum(dnbinom(o.whospi, mu=pmax(0.1, state.o.whospi$x),
                    size=whosp_nbinom_size, log=T))

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
    
    loglH.2 <- sum(dnbinom(dhospi[1:(d.reliable.hosp-1)],
                           mu=pmax(0.1, (state$y.hospi + state$o.hospi)
                                   [dstart:(dstart + d.reliable.hosp)]),
                           size=hosp_nbinom_size1, log=T)) +
        sum(dnbinom(dhospi[d.reliable.hosp:length(dhosp)],
                    mu=pmax(0.1, (state$y.hospi + state$o.hospi)
                            [(dstart + d.reliable.hosp + 1):dend]),
                    size=hosp_nbinom_size2, log=T)) +
        sum(dnbinom(dhospi[(length(dhosp) - last):length(dhosp)],
                    mu=pmax(0.1, (state$y.hospi + state$o.hospi)[(dend - last):dend]),
                    size=hosp_nbinom_size2 * f_last7, log=T))

    loglH <- loglH.1 + loglH.2
    
    ##loglH <- sum(dpois(y.whospi, lambda=pmax(0.1, state.y.whospi$x), log=T)) +
    ##         sum(dpois(o.whospi, lambda=pmax(0.1, state.o.whospi$x), log=T))

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
                           size=mort_nbinom_size * 2, log=T)) +
               sum(dnbinom(y.dmorti[(length(y.dmorti) - last):length(y.dmorti)],
                           mu=pmax(0.001, state$y.deadi[(dend - last):dend]),
                           size=mort_nbinom_size * 2 * f_last7, log=T))

    o.loglD <- sum(dnbinom(o.dmorti,
                           mu=pmax(0.001, state$o.deadi[dstart:dend]),
                           size=o.wdmorti * mort_nbinom_size, log=T)) +
               sum(dnbinom(o.dmorti[(length(o.dmorti) - last):length(o.dmorti)],
                           mu=pmax(0.001, state$o.deadi[(dend - last):dend]),
                           size=mort_nbinom_size * f_last7, log=T))
  

    result <- y.loglC + o.loglC + loglH + y.loglD + o.loglD

    if (!is.na(result)) {
        priorp <- calclogp(params)
        postp <- result + priorp
        if (postp > maxLogP) {
            maxLogP <<- postp
            print(params)
            print(c(it, y.loglC, o.loglC, loglH, y.loglD, o.loglD, priorp, result, postp))
            write.csv(c(postp, params), file="max.csv")
            state <<- calcNominalState(state)
            print(state.y.whospi$x)
            print(state.o.whospi$x)
            pdf("max.pdf", width=20, height=10)
            graphs(postp)
            dev.off()
        }
    }

    if (it %% 1000 == 0) {
        priorp <- calclogp(params)
        postp <- result + priorp
        print(params)
        print(c(it, y.loglC, o.loglC, loglH, y.loglD, o.loglD, priorp, result, postp))
        state <<- calcNominalState(state)
        graphs(postp)
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
                    "t3o", "betay3", "betao3", "betayo3",
                    "t4o", "betay4", "betao4", "betayo4",
                           "betay5", "betao5", "betayo5",
                    "t6o", "betay6", "betao6", "betayo6",
                    "t7o", "betay7", "betao7", "betayo7",
                           "betay9", "betao9", "betayo9",
                    "ifrred", "y.CL", "o.CL",
                    "t10o", "t11o", "betay11", "betao11", "betayo11")
keyparamnames <- c("betay6", "betao6", "betayo6",
                   "betay7", "betao7", "betayo7",
                   "betay9", "betao9", "betayo9",
                   "ifrred")
fitkeyparamnames <- keyparamnames

init <- c(3.1, 0.52, 0.50,
          1.5, 0.50, 0.44,
          0.8, 0.35, 0.09,
          550, 0.75, 14, 21,
          11, 1.5, 14, 21,
          14, -8, 5.2, 5.8, 0.49, 0.30,
          92, 1.0, 0.18, 0.02,
          29, 1.4, 0.4, 0.04,
              0.96, 0.27, 0.007,
          28, 1.2, 0.9, 0.02,
          32, 2.0, 0.6, 0.10,
              1.3, 0.3, 0.05,
          0.34, 6, 11,
          245, 260, 1.9, 0.3, 0.04)

df_params <- data.frame(name = fit.paramnames,
                        min = c(2 * gamma, 0.5 * gamma, 0.5 * gamma,
                                1 * gamma, 0.5 * gamma, 0.5 * gamma,
                                0.1 * gamma, 0.1 * gamma, 0.01 * gamma,
                                500, 0.5, 9, 14,
                                3, 0.5, 9, 14,
                                0, -20, 2, 2, 0.05, 0.05,
                                60, 0.5 * gamma, 0.01 * gamma, 0.002 * gamma,
                                10, 0.5 * gamma, 0.01 * gamma, 0.002 * gamma,
                                    0.1 * gamma, 0.01 * gamma, 0.002 * gamma,
                                10, 0.2 * gamma, 0.01 * gamma, 0.002 * gamma,
                                10, 0.2 * gamma, 0.01 * gamma, 0.002 * gamma,
                                    0.2 * gamma, 0.01 * gamma, 0.002 * gamma,
                                0.1, 4, 4,
                                d10 - 10, d11 - 10, 0.2 * gamma, 0.01 * gamma, 0.002 * gamma),
                        max = c(8 * gamma, 5 * gamma, 5 * gamma,
                                5 * gamma, 3 * gamma, 3 * gamma,
                                2 * gamma, 2 * gamma, 2 * gamma,
                                5000, 2, 18, 25,
                                100, 6, 18, 25,
                                max(dmort[length(dmort)] / 10, total_deaths_at_lockdown * 10),
                                10, 9, 9, 1, 1,
                                100, 3 * gamma, 1 * gamma, 0.5 * gamma,
                                50, 3 * gamma, 1 * gamma, 0.5 * gamma,
                                    2 * gamma, 1 * gamma, 0.5 * gamma,
                                50, 4 * gamma, 3 * gamma, 0.5 * gamma,
                                50, 4 * gamma, 3 * gamma, 0.5 * gamma,
                                    3 * gamma, 3 * gamma, 0.5 * gamma,
                                0.7, 14, 17,
                                duncertain - 7, duncertain, 3 * gamma, 3 * gamma, 0.5 * gamma),
                        init = init)

print(df_params)
