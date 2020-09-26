library("deSolve")
system("R CMD SHLIB model-age-simple.cpp")
dyn.load("model-age-simple.so")

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
    ydied_latency <- params[10]
    odied_latency <- params[11]
    t0_morts <- params[12]
    t0o <- params[13]
    DLsd <- params[14]
    t3o <- params[15] ## start of increase
    betay3 <- params[16]
    betao3 <- params[17]
    betayo3 <- params[18]
    t4o <- params[19] ## end of increase
    betay4 <- params[20]
    betao4 <- params[21]
    betayo4 <- params[22]
    
    ## convolution profiles
    y.died_cv_profile = diedProfile(ydied_latency, DLsd)
    o.died_cv_profile = diedProfile(odied_latency, DLsd)

    padding = max(-y.died_cv_profile$kbegin, -o.died_cv_profile$kbegin) + 1

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
               t4 = 1E10, betay4 = betay4, betao4 = betao4, betayo4 = betayo4)

    Y <- c(Sy = y.N - Initial, E1y = Initial, E2y = 0, Iy = 0, Ry = 0,
           So = o.N, E1o = 0, E2o = 0, Io = 0, Ro = 0)

    times <- (padding + 1):(padding + period)

    out <- ode(Y, times, func = "derivs", parms = parms,
               dllname = "model-age-simple",
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

        parms <- c(Ny = y.N, No = o.N,
                   a1 = a1, a2 = a2, gamma = gamma,
                   t0 = data_offset + lockdown_offset + t0o,
                   t1 = data_offset + lockdown_offset + max(t0o, 0),
                   t2 = t2,
                   betay0 = betay0, betao0 = betao0, betayo0 = betayo0,
                   betay1 = betay1, betao1 = betao1, betayo1 = betayo1,
                   betay2 = betay2, betao2 = betao2, betayo2 = betayo2,
                   t3 = t3, betay3 = betay3, betao3 = betao3, betayo3 = betayo3,
                   t4 = t4, betay4 = betay4, betao4 = betao4, betayo4 = betayo4)

        out <- ode(Y, times, func = "derivs", parms = parms,
                   dllname = "model-age-simple",
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

    posterior
}

calcNominalState <- function(state)
{
    state$died <- state$y.died + state$o.died
    state$deadi <- state$y.deadi + state$o.deadi

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
    ydied_latency <- params[10]
    odied_latency <- params[11]
    t0_morts <- params[12]
    t0o <- params[13]
    DLsd <- params[14]
    t3o <- params[15] ## start of increase
    betay3 <- params[16]
    betao3 <- params[17]
    betayo3 <- params[18]
    t4o <- params[19] ## end of increase
    betay4 <- params[20]
    betao4 <- params[21]
    betayo4 <- params[22]

    logPriorP <- 0
    
    logPriorP <- logPriorP + dnorm(t0o, mean=0, sd=10, log=T)
    logPriorP <- logPriorP + dnorm(lockdown_offset + lockdown_transition_period + t3o, mean=d3, sd=30, log=T)
    logPriorP <- logPriorP + dnorm(lockdown_offset + lockdown_transition_period + t3o + t4o, mean=d4, sd=30, log=T)
    logPriorP <- logPriorP + dnorm(betay0 - betay1, mean=0, sd=2*gamma, log=T)
    logPriorP <- logPriorP + dnorm(betao0 - betao1, mean=0, sd=2*gamma, log=T)
    logPriorP <- logPriorP + dnorm(betayo0 - betayo1, mean=0, sd=2*gamma, log=T)
    logPriorP <- logPriorP + dnorm(betay1 - betay2, mean=0, sd=2*gamma, log=T)
    logPriorP <- logPriorP + dnorm(betao1 - betao2, mean=0, sd=2*gamma, log=T)
    logPriorP <- logPriorP + dnorm(betayo1 - betayo2, mean=0, sd=2*gamma, log=T)
    logPriorP <- logPriorP + dnorm(betay3 - betay4, mean=0, sd=0.5*gamma, log=T)
    logPriorP <- logPriorP + dnorm(betao3 - betao4, mean=0, sd=0.5*gamma, log=T)
    logPriorP <- logPriorP + dnorm(betayo3 - betayo4, mean=0, sd=0.5*gamma, log=T)
    logPriorP <- logPriorP + dnorm(DLsd, mean=5, sd=1, log=T)
    logPriorP <- logPriorP + dnorm(ydied_latency, mean=21, sd=4, log=T)
    logPriorP <- logPriorP + dnorm(odied_latency, mean=21, sd=4, log=T)

    logPriorP
}

calclogl <- function(params, x) {
    state <<- calculateModel(params, FitTotalPeriod)

    if (state$offset == InvalidDataOffset)
        state$offset = 1

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

    result <- y.loglD + o.loglD
    
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
                    "y.DL",
                    "o.DL",
                    "t0_morts", "t0o",
                    "DLsd",
                    "t3o", "betay3", "betao3", "betayo3",
                    "t4o", "betay4", "betao4", "betayo4")
keyparamnames <- c("betay3", "betao3", "betayo3", "betay4", "betao4", "betayo4")
fitkeyparamnames <- keyparamnames

init <- c(3.6 * gamma, 3.6 * gamma, 3.6 * gamma,
          2.0 * gamma, 2.0 * gamma, 2.0 * gamma,
          0.8 * gamma, 0.8 * gamma, 0.8 * gamma,
          21,
          21,
          total_deaths_at_lockdown, -1,
          5,
          d3 - lockdown_offset - lockdown_transition_period, 0.8 * gamma, 0.8 * gamma, 0.8 * gamma,
          d4 - d3, 0.8 * gamma, 0.8 * gamma, 0.8 * gamma)

df_params <- data.frame(name = fit.paramnames,
                        min = c(2 * gamma, 2 * gamma, 2 * gamma,
                                1 * gamma, 1 * gamma, 1 * gamma,
                                0.05 * gamma, 0.01 * gamma, 0.01 * gamma,
                                10,
                                10,
                                0, -30,
                                2,
                                50, 0.05 * gamma, 0.01 * gamma, 0.01 * gamma,
                                60, 0.05 * gamma, 0.01 * gamma, 0.01 * gamma),
                        max = c(8 * gamma, 8 * gamma, 8 * gamma,
                                5 * gamma, 5 * gamma, 5 * gamma,
                                2 * gamma, 2 * gamma, 2 * gamma,
                                50,
                                50,
                                max(dmort[length(dmort)] / 10, total_deaths_at_lockdown * 10),
                                30,
                                9,
                                100, 3 * gamma, 3 * gamma, 3 * gamma,
                                120, 3 * gamma, 3 * gamma, 3 * gamma),
                        init = init)

print(df_params)
