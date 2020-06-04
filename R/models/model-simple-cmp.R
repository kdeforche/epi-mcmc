require(Rcpp)

InvalidDataOffset <- 10000
Initial <- 1

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

cppFunction('NumericVector odesimstepc(double S, double E, double I, double R, int N, double beta0, double a0, double Tinf0, double beta1, double a1, double Tinf1) {
   const int LOOPS = 10;
   const double Ts = 1.0/LOOPS;

   NumericVector out(6);

   for (int l = 0; l < LOOPS; ++l) {
     const double beta = 1.0/LOOPS * (l * beta1 + (LOOPS - l) * beta0);
     const double a = 1.0/LOOPS * (l * a1 + (LOOPS - l) * a0);
     const double Tinf = 1.0/LOOPS * (l * Tinf1 + (LOOPS - l) * Tinf0);

     const double gamma = 1/Tinf;

     const double inf = (beta * I) / N * S;
     const double got_infected = inf;
     const double got_infectious = a * E;
     const double got_removed = gamma * I;
   
     const double deltaS = -got_infected;
     const double deltaE = got_infected - got_infectious;
     const double deltaI = got_infectious - got_removed;
     const double deltaR = got_removed;

     if (l == 0) {
       out[4] = Tinf * inf / I;
       out[5] = out[5] * N / S;
     }

     S += deltaS * Ts;
     E += deltaE * Ts;
     I += deltaI * Ts;
     R += deltaR * Ts;
   }

   out[0] = S;
   out[1] = E;
   out[2] = I;
   out[3] = R;

   return out;
}')

simstep.C <- function(state, N, beta0, a0, Tinf0, beta1, a1, Tinf1)
{
    i = state$i

    newState <- odesimstepc(state$S[i], state$E[i], state$I[i], state$R[i], N,
                            beta0, a0, Tinf0,
                            beta1, a1, Tinf1)

    state$S[i + 1] = newState[1]
    state$E[i + 1] = newState[2]
    state$I[i + 1] = newState[3]
    state$R[i + 1] = newState[4]
    state$Re[i + 1] = newState[5]
    state$Rt[i + 1] = newState[6]

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
    Tinf0 <- params[3]
    Tinft <- params[4]
    hosp_rate <- params[5]
    hosp_latency <- params[6]
    died_latency <- params[7]
    mort_lockdown_threshold <- params[8]
    Es <- tail(params, n=-8)

    beta0 = R0 / Tinf0
    betat = Rt / Tinft
    
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
    state$hosp <- rep(0, padding + period)
    state$died <- rep(0, padding + period)
    state$Sr <- c()
    state$i <- padding + 1

    data_offset = InvalidDataOffset

    l.beta = beta0
    l.Tinf = Tinf0
    
    for (i in (padding + 1):(padding + period)) {
        beta = calcpar(i - data_offset, beta0, betat, Es)
        Tinf = calcpar(i - data_offset, Tinf0, Tinft, Es)
        
        state <- simstep.C(state, N,
                           l.beta, a, l.Tinf,
                           beta, a, Tinf)

        l.beta = beta
        l.Tinf = Tinf
        
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
    result[5] = exp(params[5])

    result
}

invTransformParams <- function(posterior)
{
    posterior$Tinc = Tinc

    posterior$beta0 = posterior$R0 / posterior$Tinf0
    posterior$betat = posterior$Rt / posterior$Tinft

    posterior$HR = exp(posterior$logHR)

    posterior
}

## log likelihood function for fitting this model to observed data:
##   dhospi and dmorti
calclogl <- function(params) {
    R0 <- params[1]
    Rt <- params[2]
    Tinf0 <- params[3]
    Tinft <- params[4]
    hosp_rate <- params[5]
    hosp_latency <- params[6]
    died_latency <- params[7]
    mort_lockdown_threshold <- params[8]
    Es <- tail(params, n=-8)

    if (R0 < 0 || Rt < 0) {
        return(-Inf)
    }
    
    if (Tinf0 < 0.2 || Tinf0 > 13.8) {
        return(-Inf)
    }

    if (Tinft < 0.2 || Tinft > 13.8) {
        return(-Inf)
    }

    beta0 = R0 / Tinf0
    G0 = Tinc + Tinf0/2
    betat = Rt / Tinft
    Gt = Tinc + Tinft/2

    if (beta0 < 0 || betat < 0) {
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
    logPriorP <- logPriorP + dnorm(died_latency, mean=10, sd=20, log=T)

    logPriorP <- logPriorP + dnorm(G0, mean=5, sd=1, log=T)
    logPriorP <- logPriorP + dnorm(G0 - Gt, mean=0, sd=1.5, log=T)

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
        print(c(R0, Rt, Tinf0, Tinft, beta0, betat))
	print(c(it, result))
	graphs()
    }

    result
}


fit.paramnames <- c("R0", "Rt", "Tinf0", "Tinft",
                    "logHR", "HL", "DL", "lockdownmort")
keyparamnames <- c("beta0", "betat", "R0", "Rt", "Tinf0", "Tinft")
fitkeyparamnames <- c("R0", "Rt", "Tinf0", "Tinft")

init <- c(2.9, 0.9, 3, 3, log(0.02), 10, 20, total_deaths_at_lockdown,
          rep(0.9, length(Es.time)))
scales <- c(1, 1, 1, 1, 0.05, 1, 1, total_deaths_at_lockdown / 20,
            rep(0.05, length(Es.time)))

if (length(Es.time) > 0) {
    for (i in 1:length(Es.time)) {
        fit.paramnames <- c(fit.paramnames, paste("E", i, sep=""))
    }
}
