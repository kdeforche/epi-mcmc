require(Rcpp)

InvalidDataOffset <- 10000
Initial <- 1

## Generation interval of 5.2 days, incubation period of 5.2, 50% presymtomatic transmission
## -> 2.5 + 2.5 days => average is 5 dagen
##    about half before 5 days, and half after ?
Tinf <- 14
Tinc <- 2.5
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

convolute <- function(values, i, profile)
{
    i1 <- i + profile$kbegin
    i2 <- i + profile$kend

    (profile$values %*% values[i1:i2])[1,1]
}

cppFunction('NumericVector odesimstepc(double S, double E, double In, double Is, double R, int N, double betaIn0, double betaIs0, double a0, double Tin0, double Tis0, double betaIn1, double betaIs1, double a1, double Tin1, double Tis1) {
   const int LOOPS = 10;
   const double Ts = 1.0/LOOPS;

   NumericVector out(7);

   for (int l = 0; l < LOOPS; ++l) {
     const double betaIn = 1.0/LOOPS * (l * betaIn1 + (LOOPS - l) * betaIn0);
     const double betaIs = 1.0/LOOPS * (l * betaIs1 + (LOOPS - l) * betaIs0);
     const double a = 1.0/LOOPS * (l * a1 + (LOOPS - l) * a0);
     const double Tin = 1.0/LOOPS * (l * Tin1 + (LOOPS - l) * Tin0);
     const double Tis = 1.0/LOOPS * (l * Tis1 + (LOOPS - l) * Tis0);

     const double gamma1 = 1/Tin;
     const double gamma2 = 1/Tis;

     const double inf1 = (betaIn * In) / N * S;
     const double inf2 = (betaIs * Is) / N * S;
     const double got_infected = inf1 + inf2;
     const double got_infectious = a * E;
     const double got_isolated = gamma1 * In;
     const double got_removed = gamma2 * Is;
   
     const double deltaS = -got_infected;
     const double deltaE = got_infected - got_infectious;
     const double deltaIn = got_infectious - got_isolated;
     const double deltaIs = got_isolated - got_removed;
     const double deltaR = got_removed;

     if (l == 0) {
       out[5] = Tin * inf1 / In + Tis * inf2 / Is;
       out[6] = out[5] * N / S;
     }

     S += deltaS * Ts;
     E += deltaE * Ts;
     In += deltaIn * Ts;
     Is += deltaIs * Ts;
     R += deltaR * Ts;
   }

   out[0] = S;
   out[1] = E;
   out[2] = In;
   out[3] = Is;
   out[4] = R;

   return out;
}')

simstep.C <- function(state, N,
                      betaIn0, betaIs0, a0, Tin0, Tis0,
                      betaIn1, betaIs1, a1, Tin1, Tis1)
{
    i = state$i

    newState <- odesimstepc(state$S[i], state$E[i], state$In[i], state$Is[i],
                            state$R[i], N,
                            betaIn0, betaIs0, a0, Tin0, Tis0,
                            betaIn1, betaIs1, a1, Tin1, Tis1)

    state$S[i + 1] = newState[1]
    state$E[i + 1] = newState[2]
    state$In[i + 1] = newState[3]
    state$Is[i + 1] = newState[4]
    state$R[i + 1] = newState[5]
    state$Re[i + 1] = newState[6]
    state$Rt[i + 1] = newState[7]

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
    R0 <- params[1]
    Rt <- params[2]
    Tef0 <- params[3]
    Teft <- params[4]
    Tin0 <- params[5]
    Tint <- params[6]
    hosp_rate <- params[7]
    hosp_latency <- params[8]
    died_latency <- params[9]
    mort_lockdown_threshold <- params[10]
    Es <- tail(params, n=-10)

    Tis0 = Tinf - Tin0
    Ris0 = R0 * (Tef0 - Tin0) / (Tinf - Tin0)
    Rin0 = R0 - Ris0
    betaIn0 = Rin0 / Tin0
    betaIs0 = Ris0 / Tis0

    Tist = Tinf - Tint
    Rist = Rt * (Teft - Tint) / (Tinf - Tint)
    Rint = Rt - Rist
    betaInt = Rint / Tint
    betaIst = Rist / Tist
    
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

    data_offset = InvalidDataOffset

    l.betaIn = betaIn0
    l.betaIs = betaIs0
    l.Tin = Tin0
    l.Tis = Tis0
    
    for (i in (padding + 1):(padding + period)) {
        betaIn = calcpar(i - data_offset, betaIn0, betaInt, Es)
        betaIs = calcpar(i - data_offset, betaIs0, betaIst, Es)
        Tin = calcpar(i - data_offset, Tin0, Tint, Es)
        Tis = Tinf - Tin
        
        state <- simstep.C(state, N,
                           l.betaIn, l.betaIs, a, l.Tin, l.Tis,
                           betaIn, betaIs, a, Tin, Tis)

        l.betaIn = betaIn
        l.betaIs = betaIs
        l.Tin = Tin
        l.Tis = Tis
        
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

    posterior$Tis0 = posterior$Tinf - posterior$Tin0
    posterior$Ris0 = posterior$R0 * (posterior$Tef0 - posterior$Tin0) /
        (posterior$Tinf - posterior$Tin0)
    posterior$Rin0 = posterior$R0 - posterior$Ris0
    posterior$betaIn0 = posterior$Rin0 / posterior$Tin0
    posterior$betaIs0 = posterior$Ris0 / posterior$Tis0
    posterior$beta0 = posterior$R0 / posterior$Tef0

    posterior$Tist = posterior$Tinf - posterior$Tint
    posterior$Rist = posterior$Rt * (posterior$Teft - posterior$Tint) /
        (posterior$Tinf - posterior$Tint)
    posterior$Rint = posterior$Rt - posterior$Rist
    posterior$betaInt = posterior$Rint / posterior$Tint
    posterior$betaIst = posterior$Rist / posterior$Tist
    posterior$betat = posterior$Rt / posterior$Teft

    posterior$deltaTeff = posterior$Teft - posterior$Tef0

    posterior$frTef = posterior$Teft / posterior$Tef0
    posterior$frbeta = posterior$betat / posterior$beta0
    posterior$frR = posterior$Rt / posterior$R0
    
    posterior$HR = exp(posterior$logHR)

    posterior
}

## log likelihood function for fitting this model to observed data:
##   dhospi and dmorti
calclogl <- function(params) {
    R0 <- params[1]
    Rt <- params[2]
    Tef0 <- params[3]
    Teft <- params[4]
    Tin0 <- params[5]
    Tint <- params[6]
    hosp_rate <- params[7]
    hosp_latency <- params[8]
    died_latency <- params[9]
    mort_lockdown_threshold <- params[10]
    Es <- tail(params, n=-10)

    if (R0 < 0 || Rt < 0) {
        return(-Inf)
    }
    
    if (Tin0 < 0.2 || Tin0 > Tinf - 0.2) {
        return(-Inf)
    }

    if (Tint < 0.2 || Tint > Tinf - 0.2) {
        return(-Inf)
    }

    if (Tef0 < 0.2 || Tef0 > Tinf - 0.2) {
        return(-Inf)
    }

    Tis0 = Tinf - Tin0
    Ris0 = R0 * (Tef0 - Tin0) / (Tinf - Tin0)
    Rin0 = R0 - Ris0
    betaIn0 = Rin0 / Tin0
    betaIs0 = Ris0 / Tis0
    beta0 = R0 / Tef0

    Tist = Tinf - Tint
    Rist = Rt * (Teft - Tint) / (Tinf - Tint)
    Rint = Rt - Rist
    betaInt = Rint / Tint
    betaIst = Rist / Tist
    betat = Rt / Teft

    if (betaIn0 < 0 || betaInt < 0 || betaIs0 < 0 || betaIst < 0) {
        return(-Inf)
    }
    
    if (betaIst > betaInt) {
        return(-Inf)
    }

    if (betaIs0 > betaIn0) {
        return(-Inf)
    }

    if (Ris0 > Rin0) {
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

    logPriorP <- logPriorP + dnorm(hosp_latency, mean=10, sd=20, log=T)
    logPriorP <- logPriorP + dnorm(died_latency, mean=18, sd=1, log=T)

    logPriorP <- logPriorP + dnorm(Tef0, mean=2.5, sd=1, log=T)
    logPriorP <- logPriorP + dnorm(Tef0 - Teft, mean=0, sd=1.5, log=T)

    ##logPriorP <- logPriorP + dnorm(betaIn0, mean=1, sd=0.1, log=T)
    ##logPriorP <- logPriorP + dnorm(betaInt, mean=1, sd=0.1, log=T)

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
        print("=========================== Increase Padding ? ===================")
        dstart <- 1
        dend <- dstart + length(dmorti) - 1
    }

    loglD <- sum(dnbinom(dmorti,
                         mu=pmax(0.1, state$deadi[dstart:dend]), size=mort_nbinom_size, log=T))

    it <<- it + 1

    result <- logPriorP + loglLD + loglH + loglD
    
    if (it %% 1000 == 0) {
        print(params)
        print(c(R0, Rt, Tef0, Teft, beta0, betat))
	print(c(it, result))
	graphs()
    }

    result
}


fit.paramnames <- c("R0", "Rt", "Tef0", "Teft", "Tin0", "Tint",
                    "logHR", "HL", "DL", "lockdownmort")
keyparamnames <- c("betaIn0", "betaInt", "betaIs0", "betaIst", "R0", "Rt", "Tef0", "Teft",
                   "Rin0", "Ris0")
fitkeyparamnames <- c("R0", "Rt", "Tef0", "Teft", "Tin0", "Tint")

init <- c(2.9, 0.9, 3, 3, 2, 1, log(0.02), 10, 20, total_deaths_at_lockdown,
          rep(0.9, length(Es.time)))
scales <- c(1, 1, 1, 1, 1, 1, 0.05, 1, 1, total_deaths_at_lockdown / 20,
            rep(0.05, length(Es.time)))

if (length(Es.time) > 0) {
    for (i in 1:length(Es.time)) {
        fit.paramnames <- c(fit.paramnames, paste("E", i, sep=""))
    }
}
