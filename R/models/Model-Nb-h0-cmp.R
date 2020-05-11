require(Rcpp)

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

cppFunction('NumericVector odesimstepc(double S, double E, double I, double R, double Ib, int N, double betaIn, double betaIs, double a, double Tinf, double h0, double BS, double prevBS) {
   double gamma = 1/Tinf;

   const int LOOPS = 10;
   const double Ts = 1.0/LOOPS;

   NumericVector out(7);

   if (prevBS != BS) {
     if (prevBS == 1) {
        Ib = h0 * Ib;
     } else {
        Ib = N / BS * (1 - pow(std::max(0.0, 1 - prevBS * Ib / N), BS/prevBS));
     }
   }

   for (int l = 0; l < LOOPS; ++l) {
     const double susceptibles_in_infected_bubbles = std::max(0.0, Ib * BS - (N - S));
     const double inf1 = (betaIn * I) / N * susceptibles_in_infected_bubbles;
     const double susceptible_bubbles = std::max(0.0, N/BS - Ib);
     const double inf2 = (betaIs * I) / N * susceptible_bubbles;
     const double got_infected = inf1 + inf2;
     const double got_infectious = a * E;
     const double got_removed = gamma * I;
   
     const double deltaS = -got_infected;
     const double deltaE = got_infected - got_infectious;
     const double deltaI = got_infectious - got_removed;
     const double deltaR = got_removed;

     if (l == 0) {
       out[5] = Tinf * got_infected / I;
       out[6] = out[5] * N / I;
     }

     S += deltaS * Ts;
     E += deltaE * Ts;
     I += deltaI * Ts;
     R += deltaR * Ts;
     Ib = std::min(0.9 * N / BS, Ib + inf2 * Ts);
   }

   out[0] = S;
   out[1] = E;
   out[2] = I;
   out[3] = R;
   out[4] = Ib;

   return out;
}')

simstep.C <- function(state, N, betaIn, betaIs, a, Tinf, h0, BS, prevBS)
{
    i = state$i

    newState <- odesimstepc(state$S[i], state$E[i], state$I[i], state$R[i],
                            state$Ib[i], N, betaIn, betaIs, a, Tinf, h0, BS, prevBS)

    state$S[i + 1] = newState[1]
    state$E[i + 1] = newState[2]
    state$I[i + 1] = newState[3]
    state$R[i + 1] = newState[4]
    state$Ib[i + 1] = newState[5]
    state$Re[i + 1] = newState[6]
    state$Rt[i + 1] = newState[7]
    state$BS[i + 1] = BS

    i = i + 1
    state$i = i

    state
}

calcIb <- function(I, N, BS, h) {
    result = (1 - (1 - h * I/N)^BS) * N / BS
    result
}

calch <- function(I, N, BS, Ib) {
    result = (1 - max(0, 1 - Ib * BS / N)^(1/BS)) * N / I
    result
}

adjustIb <- function(N, BS1, BS2, Ib) {
    result = N / BS2 * (1 - max(0, 1 - BS1 * Ib / N)^(BS2/BS1))
    result
}

simstep.R <- function(state, N, betaIn, betaIs, a, Tinf, BS, prevBS)
{
    i = state$i

    gamma = 1/Tinf

    loops = 1
    Ts = 1/loops

    S = state$S[i]
    E = state$E[i]
    I = state$I[i]
    R = state$R[i]
    Ib = state$Ib[i] ## number of infected (E+I+R) bubbles

    if (prevBS != BS) {
        ## adjust Ib keeping same heterogeneity
        Ib = adjustIb(N, prevBS, BS, Ib)
    }
    
    for (l in 1:loops) {
        susceptibles_in_bubbles = max(0, Ib * BS - (N - S))
        inf1 <- (betaIn * I) / N * susceptibles_in_bubbles
        inf2 <- (betaIs * I) / N * max(0, S - susceptibles_in_bubbles)
        got_infected = inf1 + inf2
        got_infectious = a * E
        got_removed = gamma * I

        deltaS = -got_infected
        deltaE = got_infected - got_infectious
        deltaI = got_infectious - got_removed
        deltaR = got_removed

        if (l == 1) {
            state$Re[i + 1] = Tinf * got_infected / I;
            state$Rt[i + 1] = (betaIn + betaIs) * Tinf;
        }

        S = S + deltaS * Ts
        E = E + deltaE * Ts
        I = I + deltaI * Ts
        R = R + deltaR * Ts
        Ib = min(0.9 * N / BS, Ib + inf2 * Ts)
    }

    state$S[i + 1] = S
    state$E[i + 1] = E
    state$I[i + 1] = I
    state$R[i + 1] = R
    state$Ib[i + 1] = Ib

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
    betaInt <- params[1]
    betaIs0 <- params[2]
    betaIst <- params[3]
    hosp_rate <- params[4]
    hosp_latency <- params[5]
    died_latency <- params[6]
    mort_lockdown_threshold <- params[7]
    logNbt <- params[8]
    h0 <- params[9]
    Es <- tail(params, n=-9)

    Tinf <- 8
    Tinc <- 3.5
    died_rate <- 0.007

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
    state$Ib <- rep(Initial, padding + period)
    state$BS <- rep(0, padding + period)
    state$hosp <- rep(0, padding + period)
    state$died <- rep(0, padding + period)
    state$i <- padding + 1

    data_offset = InvalidDataOffset

    if (TRUE) {
        prevBS = 1
        for (i in (padding + 1):(padding + period)) {
            betaIn = betaInt ## calcpar(i - data_offset, 0, betaInt, Es)
            betaIs = calcpar(i - data_offset, betaIs0, betaIst, Es)
            logNb = calcpar(i - data_offset, log(N), logNbt, Es)

            if (logNb == log(N))
                Nb = N
            else
                Nb = exp(logNb)
            
            BS = N / Nb
            
            state <- simstep.C(state, N, betaIn, betaIs, a, Tinf, h0, BS, prevBS)
            prevBS = BS

            s = convolute(state$S, i, hosp_cv_profile)
            state$hosp[i] <- (N - s) * hosp_rate

            ## died:
            ##  from those that (ever) became infectious           
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
    result[4] = exp(params[4])

    result
}

invTransformParams <- function(posterior)
{
    posterior$HR = exp(posterior$logHR)
    posterior$Tinf = 8
    posterior$Nbt = exp(posterior$logNbt)

    posterior$R0 = posterior$betaIs0 * posterior$Tinf

    posterior
}

## log likelihood function for fitting this model to observed data:
##   dhospi and dmorti
calclogl <- function(params) {
    betaInt <- params[1]
    betaIs0 <- params[2]
    betaIst <- params[3]
    hosp_rate <- params[4]
    hosp_latency <- params[5]
    died_latency <- params[6]
    mort_lockdown_threshold <- params[7]
    logNbt <- params[8]
    h0 <- params[9]
    Es <- tail(params, n=-9)

    if (betaInt < 1E-10 || betaIs0 < 1E-10 || betaIst < 1E-10) {
        return(-Inf)
    }

    if (logNbt < 0 || logNbt > log(N)) {
        return(-Inf)
    }

    if (h0 < 0 || h0 > 1) {
        return (-Inf)
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

    ##logPriorP <- logPriorP + dnorm(logNbt, mean=log(N), sd=4, log=T)
    ##logPriorP <- logPriorP + dnorm(betaInt, mean=0, sd=1, log=T)
    ##logPriorP <- logPriorP + dnorm(betaIs0, mean=0, sd=1, log=T)
    ##logPriorP <- logPriorP + dnorm(betaIst, mean=0, sd=1, log=T)
    logPriorP <- logPriorP + dnorm(hosp_latency, mean=10, sd=20, log=T)
    logPriorP <- logPriorP + dnorm(died_latency, mean=10, sd=20, log=T)

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

    di <- dhospi

    if (dstart < 1) {
        print("=========================== Increase Padding ? ===================")
        dstart <- 1
        dend <- dstart + length(dhospi) - 1
    }
    
    loglH <- sum(dnbinom(dhospi,
                         mu=pmax(0.1, state$hospi[dstart:dend]),
                         size=hosp_nbinom_size, log=T))

    dstart <- state$offset
    dend <- state$offset + length(dmorti) - 1

    if (dend > length(state$deadi)) {
        print("=========================== Increase FitTotalPeriod ===================")
        dend <- length(state$deadi)
        dstart <- dend - length(dmorti) + 1
    }

    if (dstart < 1) {
        print(paste(state$offset, "=========================== Increase Padding ? ==================="))
        dstart <- 1
        dend <- dstart + length(dmorti) - 1
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

fit.paramnames <- c("betaInt", "betaIs0", "betaIst", "logHR", "HL", "DL",
                    "lockdownmort", "logNbt", "h0")
keyparamnames <- c("betaInt", "betaIs0", "betaIst", "R0", "Rt", "Nbt", "h0")
fitkeyparamnames <- c("betaInt", "betaIs0", "betaIst", "logNbt", "h0")

init <- c(1, 1, 0.1, log(0.05), 15, 20, total_deaths_at_lockdown, log(N) - 1, 0.7,
          rep(0.9, length(Es.time)))
scales <- c(0.3, 0.3, 0.3, 0.05, 1, 1, total_deaths_at_lockdown / 20, 1, 0.1,
            rep(0.05, length(Es.time)))

if (length(Es.time) > 0) {
    for (i in 1:length(Es.time)) {
        fit.paramnames <- c(fit.paramnames, paste("E", i, sep=""))
    }
}
