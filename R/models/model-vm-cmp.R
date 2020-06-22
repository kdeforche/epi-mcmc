require(Rcpp)

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

     const double Tinf = 8.0;
     const double gamma = 1/Tinf;
     const double gamman = 1/Tin;

     const double inf1 = (betaIn * In) / N * S;
     const double inf2 = (betaIs * Is) / N * S;
     const double got_infected = inf1 + inf2;
     const double got_infectious = a * E;
     const double got_isolated = gamman * In;
     const double got_removed = gamma * Is;
     const double didnt_isolate = gamma * In;
   
     const double deltaS = -got_infected;
     const double deltaE = got_infected - got_infectious;
     const double deltaIn = got_infectious - got_isolated - didnt_isolate;
     const double deltaIs = got_isolated - got_removed;
     const double deltaR = got_removed + didnt_isolate;

     if (l == 0) {
       out[5] = 1 / (gamma + gamman) * inf1 / In + Tinf * inf2 / Is;
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
    pari = time - lockdown_offset + 1

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
    G0 <- params[3]
    Gt <- params[4]
    Tin0 <- params[5]
    Tint <- params[6]
    hosp_rate <- params[7]
    hosp_latency <- params[8]
    died_latency <- params[9]
    mort_lockdown_threshold <- params[10]
    Es <- tail(params, n=-10)

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

    data_offset = InvalidDataOffset

    l.betaIn = betaIn0
    l.betaIs = betaIs0
    l.Tin = Tin0
    
    for (i in (padding + 1):(padding + period)) {
        betaIn = calcpar(i - data_offset, betaIn0, betaInt, Es)
        betaIs = calcpar(i - data_offset, betaIs0, betaIst, Es)
        Tin = calcpar(i - data_offset, Tin0, Tint, Es)
        
        state <- simstep.C(state, N,
                           l.betaIn, l.betaIs, a, l.Tin, 0,
                           betaIn, betaIs, a, Tin, 0)

        l.betaIn = betaIn
        l.betaIs = betaIs
        l.Tin = Tin
        
        s1 = convolute(state$S, i, hosp_cv_profile)
        state$hosp[i] <- (N - s1) * hosp_rate

        s2 = convolute(state$S, i, died_cv_profile)
        state$died[i] <- (N - s2) * died_rate

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

## log likelihood function for fitting this model to observed data:
##   dhospi and dmorti
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

    if (G0 - Tinc < 0.2) {
        return(-Inf)
    }

    if (Gt - Tinc < 0.2) {
        return(-Inf)
    }

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
    logPriorP <- logPriorP + dnorm(died_latency, mean=21, sd=5, log=T)

    logPriorP <- logPriorP + dnorm(G0, mean=5, sd=1, log=T)
    logPriorP <- logPriorP + dnorm(G0 - Gt, mean=0, sd=1.5, log=T)

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
        print(c(R0, Rt, G0, Gt, beta0, betat))
	print(c(it, result))
	graphs()
    }

    result
}


fit.paramnames <- c("R0", "Rt", "G0", "Gt", "Tin0", "Tint",
                    "logHR", "HL", "DL", "lockdownmort")
keyparamnames <- c("betaIn0", "betaInt", "betaIs0", "betaIst", "R0", "Rt", "G0", "Gt",
                   "Rin0", "Ris0")
fitkeyparamnames <- c("R0", "Rt", "G0", "Gt", "Tin0", "Tint")

init <- c(2.9, 0.9, 5, 5, 2, 1, log(0.02), 10, 20, total_deaths_at_lockdown,
          rep(0.9, length(Es.time)))
scales <- c(1, 1, 1, 1, 1, 1, 0.05, 1, 1, total_deaths_at_lockdown / 20,
            rep(0.05, length(Es.time)))

if (length(Es.time) > 0) {
    for (i in 1:length(Es.time)) {
        fit.paramnames <- c(fit.paramnames, paste("E", i, sep=""))
    }
}
