library("deSolve")
system("R CMD SHLIB model-age3.cpp")
dyn.load("model-age3.so")

InvalidDataOffset <- 10000
Initial <- 1

G <- 5.2
Tinc <- 3
y.died_rate <- 0.00014 / 2
m.died_rate <- 0.0028 / 2
o.died_rate <- 0.04 / 2
a <- 1/Tinc
gamma <- 1/((G - Tinc) * 2)

HLsd <- 4.1
DLsd <- 6.5

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

convolute <- function(values, i1, i2, profile)
{
    filter(values, profile$values, method="convolution", sides=1)[(i1 + profile$kend):(i2 + profile$kend)]
}

hospProfile <- calcGammaProfile
diedProfile <- calcGammaProfile

calculateModel <- function(params, period)
{
    yhosp_rate <- params[1]
    yhosp_latency <- params[2]
    ydied_latency <- params[3]
    mhosp_rate <- params[4]
    mhosp_latency <- params[5]
    mdied_latency <- params[6]
    ohosp_rate <- params[7]
    ohosp_latency <- params[8]
    odied_latency <- params[9]
    t0_morts <- params[10]

    betay0 <- params[11]
    betam0 <- params[12]
    betao0 <- params[13]
    betaym0 <- params[14]
    betayo0 <- params[15]
    betamy0 <- params[16]
    betamo0 <- params[17]

    betay1 <- params[18]
    betam1 <- params[19]
    betao1 <- params[20]
    betaym1 <- params[21]
    betayo1 <- params[22]
    betamy1 <- params[23]
    betamo1 <- params[24]

    t2o <- params[25]
    betay2 <- params[26]
    betam2 <- params[27]
    betao2 <- params[28]
    betaym2 <- params[29]
    betayo2 <- params[30]
    betamy2 <- params[31]
    betamo2 <- params[32]

    t3o <- params[33]
    betay3 <- params[34]
    betam3 <- params[35]
    betao3 <- params[36]
    betaym3 <- params[37]
    betayo3 <- params[38]
    betamy3 <- params[39]
    betamo3 <- params[40]

    ## t4o <- params[41]
    ## betay4 <- params[42]
    ## betam4 <- params[43]
    ## betao4 <- params[44]
    ## betaym4 <- params[45]
    ## betayo4 <- params[46]
    ## betamy4 <- params[47]
    ## betamo4 <- params[48]

    ## convolution profile to infer hospitalisation count
    y.hosp_cv_profile = hospProfile(yhosp_latency, HLsd)
    m.hosp_cv_profile = hospProfile(mhosp_latency, HLsd)
    o.hosp_cv_profile = hospProfile(ohosp_latency, HLsd)
    y.died_cv_profile = diedProfile(ydied_latency, DLsd)
    m.died_cv_profile = diedProfile(mdied_latency, DLsd)
    o.died_cv_profile = diedProfile(odied_latency, DLsd)

    padding = max(-y.hosp_cv_profile$kbegin, -y.died_cv_profile$kbegin,
                  -m.hosp_cv_profile$kbegin, -m.died_cv_profile$kbegin,
                  -o.hosp_cv_profile$kbegin, -o.died_cv_profile$kbegin) + 1

    state <- NULL

    state$y.S <- rep(y.N - Initial, padding + period)
    state$y.E <- rep(Initial, padding + period)
    state$y.I <- rep(0, padding + period)
    state$y.R <- rep(0, padding + period)
    state$y.hosp <- rep(0, padding + period)
    state$y.died <- rep(0, padding + period)

    state$m.S <- rep(m.N - Initial, padding + period)
    state$m.E <- rep(Initial, padding + period)
    state$m.I <- rep(0, padding + period)
    state$m.R <- rep(0, padding + period)
    state$m.hosp <- rep(0, padding + period)
    state$m.died <- rep(0, padding + period)

    state$o.S <- rep(o.N, padding + period)
    state$o.E <- rep(0, padding + period)
    state$o.I <- rep(0, padding + period)
    state$o.R <- rep(0, padding + period)
    state$o.hosp <- rep(0, padding + period)
    state$o.died <- rep(0, padding + period)

    state$Re <- rep(0, padding + period)
    state$Rt <- rep(0, padding + period)
    state$y.Re <- rep(0, padding + period)
    state$m.Re <- rep(0, padding + period)
    state$o.Re <- rep(0, padding + period)
    
    state$i <- padding + 1

    betay4 = betay3
    betam4 = betam3
    betao4 = betao3
    betaym4 = betaym3
    betayo4 = betayo3
    betamy4 = betamy3
    betamo4 = betamo3
    
    parms <- c(Ny = y.N, Nm = m.N, No = o.N,
               a = a, gamma = gamma,
               t0 = 1E10, t1 = 1E10,
               betay0 = betay0, betam0 = betam0, betao0 = betao0,
               betaym0 = betaym0, betayo0 = betayo0, betamy0 = betamy0, betamo0 = betamo0,
               betay1 = betay1, betam1 = betam1, betao1 = betao1,
               betaym1 = betaym1, betayo1 = betayo1, betamy1 = betamy1, betamo1 = betamo1,
               t2 = 1E10, betay2 = betay2, betam2 = betam2, betao2 = betao2,
               betaym2 = betaym2, betayo2 = betayo2, betamy2 = betamy2, betamo2 = betamo2,
               t3 = 1E10, betay3 = betay3, betam3 = betam3, betao3 = betao3,
               betaym3 = betaym3, betayo3 = betayo3, betamy3 = betamy3, betamo3 = betamo3,
               t4 = 1E10, betay4 = betay4, betam4 = betam4, betao4 = betao4,
               betaym4 = betaym4, betayo4 = betayo4, betamy4 = betamy4, betamo4 = betamo4)

    Y <- c(Sy = y.N - Initial, Ey = Initial, Iy = 0, Ry = 0,
           Sm = m.N - Initial, Em = Initial, Im = 0, Rm = 0,
           So = o.N, Eo = 0, Io = 0, Ro = 0)

    times <- (padding + 1):(padding + period)

    out <- ode(Y, times, func = "derivs", parms = parms,
               dllname = "model-age3",
               initfunc = "initmod", nout = 5, outnames = c("Re", "Rt", "y.Re", "m.Re", "o.Re"))

    state$y.S[(padding + 1):(padding + period)] = out[,2]
    state$m.S[(padding + 1):(padding + period)] = out[,6]
    state$o.S[(padding + 1):(padding + period)] = out[,10]

    s2 <- convolute(state$y.S, padding + 1, padding + period, y.died_cv_profile)
    state$y.died[(padding + 1):(padding + period)] = (y.N - s2) * y.died_rate

    s2 <- convolute(state$m.S, padding + 1, padding + period, m.died_cv_profile)
    state$m.died[(padding + 1):(padding + period)] = (m.N - s2) * m.died_rate

    s2 <- convolute(state$o.S, padding + 1, padding + period, o.died_cv_profile)
    state$o.died[(padding + 1):(padding + period)] = (o.N - s2) * o.died_rate

    state$died = state$y.died + state$m.died + state$o.died

    data_offset = InvalidDataOffset

    lds <- which(state$died > t0_morts)
    if (length(lds) > 0) {
        data_offset <- lds[1] - lockdown_offset

        t0 <- data_offset + lockdown_offset
        t1 <- data_offset + lockdown_offset + lockdown_transition_period
        t2 <- t1 + t2o
        t3 <- t2 + t3o
        t4 <- 1E10
        
        parms <- c(Ny = y.N, Nm = m.N, No = o.N,
                   a = a, gamma = gamma,
                   t0 = t0, t1 = t1,
                   betay0 = betay0, betam0 = betam0, betao0 = betao0,
                   betaym0 = betaym0, betayo0 = betayo0, betamy0 = betamy0, betamo0 = betamo0,
                   betay1 = betay1, betam1 = betam1, betao1 = betao1,
                   betaym1 = betaym1, betayo1 = betayo1, betamy1 = betamy1, betamo1 = betamo1,
                   t2 = t2, betay2 = betay2, betam2 = betam2, betao2 = betao2,
                   betaym2 = betaym2, betayo2 = betayo2, betamy2 = betamy2, betamo2 = betamo2,
                   t3 = t3, betay3 = betay3, betam3 = betam3, betao3 = betao3,
                   betaym3 = betaym3, betayo3 = betayo3, betamy3 = betamy3, betamo3 = betamo3,
                   t4 = t4, betay4 = betay4, betam4 = betam4, betao4 = betao4,
                   betaym4 = betaym4, betayo4 = betayo4, betamy4 = betamy4, betamo4 = betamo4)

        out <- ode(Y, times, func = "derivs", parms = parms,
                   dllname = "model-age3",
                   initfunc = "initmod", nout = 5, outnames = c("Re", "Rt", "y.Re", "m.Re", "o.Re"))
     }

    state$y.S[(padding + 1):(padding + period)] = out[,2]
    state$y.E[(padding + 1):(padding + period)] = out[,3]
    state$y.I[(padding + 1):(padding + period)] = out[,4]
    state$y.R[(padding + 1):(padding + period)] = out[,5]
    state$m.S[(padding + 1):(padding + period)] = out[,6]
    state$m.E[(padding + 1):(padding + period)] = out[,7]
    state$m.I[(padding + 1):(padding + period)] = out[,8]
    state$m.R[(padding + 1):(padding + period)] = out[,9]
    state$o.S[(padding + 1):(padding + period)] = out[,10]
    state$o.E[(padding + 1):(padding + period)] = out[,11]
    state$o.I[(padding + 1):(padding + period)] = out[,12]
    state$o.R[(padding + 1):(padding + period)] = out[,13]
    state$Re[(padding + 1):(padding + period)] = out[,14]
    state$Rt[(padding + 1):(padding + period)] = out[,15]
    state$y.Re[(padding + 1):(padding + period)] = out[,16]
    state$m.Re[(padding + 1):(padding + period)] = out[,17]
    state$o.Re[(padding + 1):(padding + period)] = out[,18]

    s1 <- convolute(state$y.S, padding + 1, padding + period, y.hosp_cv_profile)
    state$y.hosp[(padding + 1):(padding + period)] = (y.N - s1) * yhosp_rate
    
    s2 <- convolute(state$y.S, padding + 1, padding + period, y.died_cv_profile)
    state$y.died[(padding + 1):(padding + period)] = (y.N - s2) * y.died_rate

    state$y.deadi <- c(state$y.died[1], diff(state$y.died))
    state$y.hospi <- c(state$y.hosp[1], diff(state$y.hosp))

    s1 <- convolute(state$m.S, padding + 1, padding + period, m.hosp_cv_profile)
    state$m.hosp[(padding + 1):(padding + period)] = (m.N - s1) * mhosp_rate
    
    s2 <- convolute(state$m.S, padding + 1, padding + period, m.died_cv_profile)
    state$m.died[(padding + 1):(padding + period)] = (m.N - s2) * m.died_rate

    state$m.deadi <- c(state$m.died[1], diff(state$m.died))
    state$m.hospi <- c(state$m.hosp[1], diff(state$m.hosp))

    s1 <- convolute(state$o.S, padding + 1, padding + period, o.hosp_cv_profile)
    state$o.hosp[(padding + 1):(padding + period)] = (o.N - s1) * ohosp_rate
    
    s2 <- convolute(state$o.S, padding + 1, padding + period, o.died_cv_profile)
    state$o.died[(padding + 1):(padding + period)] = (o.N - s2) * o.died_rate

    state$o.deadi <- c(state$o.died[1], diff(state$o.died))
    state$o.hospi <- c(state$o.hosp[1], diff(state$o.hosp))

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

    posterior
}

calcNominalState <- function(state)
{
    state$hosp <- state$y.hosp + state$m.hosp + state$o.hosp
    state$died <- state$y.died + state$m.died + state$o.died
    state$hospi <- state$y.hospi + state$m.hospi + state$o.hospi
    state$deadi <- state$y.deadi + state$m.deadi + state$o.deadi

    state
}

## log likelihood function for fitting this model to observed data:
##   y.dhospi, o.dhospi, y.dmorti, o.dmorti
calclogp <- function(params) {
    yhosp_rate <- params[1]
    yhosp_latency <- params[2]
    ydied_latency <- params[3]
    mhosp_rate <- params[4]
    mhosp_latency <- params[5]
    mdied_latency <- params[6]
    ohosp_rate <- params[7]
    ohosp_latency <- params[8]
    odied_latency <- params[9]
    t0_morts <- params[10]

    betay0 <- params[11]
    betam0 <- params[12]
    betao0 <- params[13]
    betaym0 <- params[14]
    betayo0 <- params[15]
    betamy0 <- params[16]
    betamo0 <- params[17]

    betay1 <- params[18]
    betam1 <- params[19]
    betao1 <- params[20]
    betaym1 <- params[21]
    betayo1 <- params[22]
    betamy1 <- params[23]
    betamo1 <- params[24]

    t2o <- params[25]
    betay2 <- params[26]
    betam2 <- params[27]
    betao2 <- params[28]
    betaym2 <- params[29]
    betayo2 <- params[30]
    betamy2 <- params[31]
    betamo2 <- params[32]

    t3o <- params[33]
    betay3 <- params[34]
    betam3 <- params[35]
    betao3 <- params[36]
    betaym3 <- params[37]
    betayo3 <- params[38]
    betamy3 <- params[39]
    betamo3 <- params[40]

    ## t4o <- params[41]
    ## betay4 <- params[42]
    ## betam4 <- params[43]
    ## betao4 <- params[44]
    ## betaym4 <- params[45]
    ## betayo4 <- params[46]
    ## betamy4 <- params[47]
    ## betamo4 <- params[48]

    logPriorP <- 0
    
    logPriorP <- logPriorP + dnorm(ydied_latency, mean=21, sd=4, log=T)
    logPriorP <- logPriorP + dnorm(mdied_latency, mean=21, sd=4, log=T)
    logPriorP <- logPriorP + dnorm(odied_latency, mean=21, sd=4, log=T)

    logPriorP <- logPriorP + dnorm(t2o, mean=(d5 - lockdown_offset - lockdown_transition_period - 15),
                                   sd=10, log=T)
    logPriorP <- logPriorP + dnorm(t3o, mean=15, sd=10, log=T)

    logPriorP <- logPriorP + dnorm(betay0, mean=0.5, sd=0.4, log=T)    
    logPriorP <- logPriorP + dnorm(betam0, mean=0.5, sd=0.4, log=T)    
    logPriorP <- logPriorP + dnorm(betao0, mean=0.5, sd=0.4, log=T)    
    logPriorP <- logPriorP + dnorm(betaym0, mean=0.5, sd=0.4, log=T)
    logPriorP <- logPriorP + dnorm(betayo0, mean=0.5, sd=0.4, log=T)
    logPriorP <- logPriorP + dnorm(betamy0, mean=0.5, sd=0.4, log=T)
    logPriorP <- logPriorP + dnorm(betamo0, mean=0.5, sd=0.4, log=T)

    logPriorP <- logPriorP + dnorm(betay1, mean=0.1, sd=0.05, log=T)    
    logPriorP <- logPriorP + dnorm(betam1, mean=0.1, sd=0.05, log=T)    
    logPriorP <- logPriorP + dnorm(betao1, mean=0.1, sd=0.05, log=T)    
    logPriorP <- logPriorP + dnorm(betaym1, mean=0.1, sd=0.05, log=T)
    logPriorP <- logPriorP + dnorm(betayo1, mean=0.1, sd=0.05, log=T)
    logPriorP <- logPriorP + dnorm(betamy1, mean=0.1, sd=0.05, log=T)
    logPriorP <- logPriorP + dnorm(betamo1, mean=0.1, sd=0.05, log=T)

    logPriorP <- logPriorP + dnorm(betay1 - betay2, mean=0, sd=0.1, log=T)
    logPriorP <- logPriorP + dnorm(betam1 - betam2, mean=0, sd=0.1, log=T)
    logPriorP <- logPriorP + dnorm(betao1 - betao2, mean=0, sd=0.1, log=T)
    logPriorP <- logPriorP + dnorm(betaym1 - betaym2, mean=0, sd=0.1, log=T)
    logPriorP <- logPriorP + dnorm(betayo1 - betayo2, mean=0, sd=0.1, log=T)
    logPriorP <- logPriorP + dnorm(betamy1 - betamy2, mean=0, sd=0.1, log=T)
    logPriorP <- logPriorP + dnorm(betamo1 - betamo2, mean=0, sd=0.1, log=T)

    logPriorP <- logPriorP + dnorm(betay2 - betay3, mean=0, sd=0.1, log=T)
    logPriorP <- logPriorP + dnorm(betam2 - betam3, mean=0, sd=0.1, log=T)
    logPriorP <- logPriorP + dnorm(betao2 - betao3, mean=0, sd=0.1, log=T)
    logPriorP <- logPriorP + dnorm(betaym2 - betaym3, mean=0, sd=0.1, log=T)
    logPriorP <- logPriorP + dnorm(betayo2 - betayo3, mean=0, sd=0.1, log=T)
    logPriorP <- logPriorP + dnorm(betamy2 - betamy3, mean=0, sd=0.1, log=T)
    logPriorP <- logPriorP + dnorm(betamo2 - betamo3, mean=0, sd=0.1, log=T)

    logPriorP
}

calclogl <- function(params, x) {
    state <<- calculateModel(params, FitTotalPeriod)

    if (state$offset == InvalidDataOffset)
        state$offset = 1

    dstart <- state$offset
    dend <- state$offset + length(y.dhospi) - 1

    if (dend > length(state$y.hospi)) {
        ##print("=========================== Increase FitTotalPeriod ===================")
        dend <- length(state$y.hospi)
        dstart <- dend - length(y.dhospi) + 1
    }

    di <- y.dhospi

    if (dstart < 1) {
        ##print("=========================== Increase Padding ? ===================")
        dstart <- 1
        dend <- dstart + length(y.dhospi) - 1
    }

    y.loglH <- sum(dnbinom(y.dhospi[1:(d.reliable.cases-1)],
                           mu=pmax(0.1, state$y.hospi[dstart:(dstart + d.reliable.cases)]),
                           size=hosp_nbinom_size1, log=T)) +
               sum(dnbinom(y.dhospi[d.reliable.cases:length(y.dhospi)],
                           mu=pmax(0.1, state$y.hospi[(dstart + d.reliable.cases + 1):dend]),
                           size=hosp_nbinom_size2, log=T))

    m.loglH <- sum(dnbinom(m.dhospi[1:(d.reliable.cases-1)],
                           mu=pmax(0.1, state$m.hospi[dstart:(dstart + d.reliable.cases)]),
                           size=hosp_nbinom_size1, log=T)) +
               sum(dnbinom(m.dhospi[d.reliable.cases:length(m.dhospi)],
                           mu=pmax(0.1, state$m.hospi[(dstart + d.reliable.cases + 1):dend]),
                           size=hosp_nbinom_size2, log=T))

    o.loglH <- sum(dnbinom(o.dhospi[1:(d.reliable.cases-1)],
                           mu=pmax(0.1, state$o.hospi[dstart:(dstart + d.reliable.cases)]),
                           size=hosp_nbinom_size1, log=T)) +
               sum(dnbinom(o.dhospi[d.reliable.cases:length(o.dhospi)],
                           mu=pmax(0.1, state$o.hospi[(dstart + d.reliable.cases + 1):dend]),
                           size=hosp_nbinom_size2, log=T))

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
                           mu=pmax(0.1, state$y.deadi[dstart:dend]),
                           size=mort_nbinom_size, log=T))

    m.loglD <- sum(dnbinom(m.dmorti,
                           mu=pmax(0.1, state$m.deadi[dstart:dend]),
                           size=mort_nbinom_size, log=T))

    o.loglD <- sum(dnbinom(o.dmorti,
                           mu=pmax(0.1, state$o.deadi[dstart:dend]),
                           size=mort_nbinom_size, log=T))

    it <<- it + 1

    result <- y.loglH + m.loglH + o.loglH + y.loglD + m.loglD + o.loglD
    
    if (it %% 1000 == 0) {
        print(params)
	print(c(it, result))
        state <<- calcNominalState(state)
	graphs()
    }

    result
}

fit.paramnames <- c("yhosp_rate", "yhosp_latency",
                    "ydied_latency", "mhosp_rate",
                    "mhosp_latency", "mdied_latency",
                    "ohosp_rate", "ohosp_latency",
                    "odied_latency", "t0_morts",
                    "betay0", "betam0", "betao0",
                    "betaym0", "betayo0", "betamy0", "betamo0",
                    "betay1", "betam1", "betao1",
                    "betaym1", "betayo1", "betamy1", "betamo1",
                    "t2o",
                    "betay2", "betam2", "betao2",
                    "betaym2", "betayo2", "betamy2", "betamo2",
                    "t3o",
                    "betay3", "betam3", "betao3",
                    "betaym3", "betayo3", "betamy3", "betamo3")
keyparamnames <- c("betay3", "betam3", "betao3", "betaym3", "betayo3")
fitkeyparamnames <- keyparamnames

init <- c(0.11, 10, 21,
          0.11, 10, 21,
          0.2, 10, 21,
          10,
          0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
          0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
          80, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
          15, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2)

df_params <- data.frame(name = fit.paramnames,
                        min = c(0.1, 7, 16,
                                0.1, 7, 16,
                                0.1, 7, 16,
                                0,
                                0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                                0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015,
                                60, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015,
                                10, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015),
                        max = c(0.3, 12, 24,
                                0.3, 12, 24,
                                0.4, 12, 24,
                                30,
                                1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5,
                                0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                                110, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                                40, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4),
                        init = init)
