require(Rcpp)
require(ggplot2)
require(gridExtra)

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

cppFunction('NumericVector odesimstepc(double S, double E, double In, double Is, double R, int N, double betaIn0, double betaIs0, double a0, double Tin0, double Tis0, double betaIn1, double betaIs1, double a1, double Tin1, double Tis1) {
   const int LOOPS = 1;
   const double Ts = 1.0/10;

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
    
simstep.R <- function(state, N, betaIn, betaIs, a, Tin, Tis, k)
{
    i = state$i

    gamma1 = 1/Tin
    gamma2 = 1/Tis

    loops = 10
    Ts = 1/loops

    S = state$S[i]
    E = state$E[i]
    In = state$In[i]
    Is = state$Is[i]
    R = state$R[i]

    day_infected1 <- 0
    day_infected2 <- 0

    for (l in 1:loops) {
        EIR = N - S
        Ne = max(EIR, k * N)

        inf1 <- (betaIn * In) / N * (Ne - EIR)
        inf2 <- (betaIs * Is) / N * S
        got_infected = inf1 + inf2
        got_infectious = a * E
        got_isolated = gamma1 * In
        got_removed = gamma2 * Is

        ##print(c(i, Tin, Tis, got_infected, got_infectious, got_isolated, got_removed))
    
        deltaS = -got_infected
        deltaE = got_infected - got_infectious
        deltaIn = got_infectious - got_isolated
        deltaIs = got_isolated - got_removed
        deltaR = got_removed

        S = S + deltaS * Ts
        E = E + deltaE * Ts
        In = In + deltaIn * Ts
        Is = Is + deltaIs * Ts
        R = R + deltaR * Ts

        if (l == 1) {
            day_infected1 <- inf1
            day_infected2 <- inf2
        }
    }

    state$Re[i + 1] = (day_infected1 + day_infected2) / ((state$In[i] + state$Is[i]) / (Tin + Tis))

    state$Rt[i + 1] = betaIn * Tin + betaIs * Tis

    state$S[i + 1] = S
    state$E[i + 1] = E
    state$In[i + 1] = In
    state$Is[i + 1] = Is
    state$R[i + 1] = R

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
    betaIn0 <- params[1]
    betaInt <- params[2]
    betaIs0 <- params[3]
    betaIst <- params[4]
    hosp_rate <- params[5]
    hosp_latency <- params[6]
    died_latency <- params[7]
    mort_lockdown_threshold <- params[8]
    gamma.in0 <- params[9]
    gamma.int <- params[10]
    Es <- tail(params, n=-11)

    Tinf <- 8
    Tinc <- 2.4
    died_rate <- 0.007

    h0 <- 1
    a <- 1 / Tinc

    os = 10
    
    hosp_cv_profile = calcConvolveProfile(-hosp_latency * os, 5 * os)
    died_cv_profile = calcConvolveProfile(-died_latency * os, 5 * os)

    padding = max(-hosp_cv_profile$kbegin, -died_cv_profile$kbegin) + 1
    total = padding + os*period
    state <- NULL
    state$Re <- rep(0, total)
    state$Rt <- rep(0, total)
    state$S <- rep(N-Initial, total)
    state$E <- rep(Initial, total)
    state$In <- rep(0, total)
    state$Is <- rep(0, total)
    state$R <- rep(0, total)
    state$hosp <- rep(0, total)
    state$died <- rep(0, total)
    state$Sr <- c()

    data_offset = InvalidDataOffset

    l.betaIn = betaIn0
    l.betaIs = betaIs0
    l.Tin = 1/gamma.in0
    l.Tis = Tinf - l.Tin

    data_offset <- 40
    
    for (step in (padding + 1):(padding + os*period)) {
        state$i <- step
        i = step / os
        betaIn = calcpar(i - data_offset, betaIn0, betaInt, Es)
        betaIs = calcpar(i - data_offset, betaIs0, betaIst, Es)
        Tin = calcpar(i - data_offset, 1/gamma.in0, 1/gamma.int, Es)
        Tis = Tinf - Tin

        state <- simstep.C(state, N,
                           l.betaIn, l.betaIs, a, l.Tin, l.Tis,
                           betaIn, betaIs, a, Tin, Tis)

        l.betaIn = betaIn
        l.betaIs = betaIs
        l.Tin = Tin
        l.Tis = Tis
        
        s = convolute(state$S, step, hosp_cv_profile)
        state$hosp[step] <- (N - s) * hosp_rate

        ## died:
        ##  from those that (ever) became infectious           
        r = convolute(state$In + state$Is + state$R, step, died_cv_profile)
        state$died[step] <- r * died_rate

        ## assuming lockdown decision was based on a cumulative mort count, with some
        ## uncertainty on the exact value due to observed cumulative mort count being a
        ## sample
        if (data_offset == InvalidDataOffset && state$died[step] > mort_lockdown_threshold) {
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
    result[1] = params[1] * result[9]
    result[2] = params[2] * result[10]
    result[3] = (params[3] - params[1]) / (8 - 1 / result[9])
    result[4] = (params[4] - params[2]) / (8 - 1 / result[10])
    result[5] = exp(params[5])

    result
}

invTransformParams <- function(posterior)
{
    posterior$RIs0 = posterior$R0 - posterior$RIn0
    posterior$RIst = posterior$Rt - posterior$RInt
    posterior$betaIn0 = posterior$RIn0 * posterior$gamma.in0 
    posterior$betaInt = posterior$RInt * posterior$gamma.int
    posterior$Tin0 = 1/posterior$gamma.in0
    posterior$Tint = 1/posterior$gamma.int
    posterior$betaIs0 = posterior$RIs0 / (8 - posterior$Tin0)
    posterior$betaIst = posterior$RIst / (8 - posterior$Tint)
    posterior$Tef0 = (posterior$Tin0 * posterior$RIn0 + 8 * posterior$RIs0) / posterior$R0
    posterior$Teft = (posterior$Tint * posterior$RInt + 8 * posterior$RIst) / posterior$Rt
    posterior$beta0 = posterior$R0 / posterior$Tef0
    posterior$betat = posterior$Rt / posterior$Teft
    ##    (posterior$Tin0 * posterior$betaIn0 + (8 - posterior$Tin0) * posterior$betaIs0)
    
    posterior$HR = exp(posterior$logHR)
    posterior$Tinf = 8

    ##posterior$R0 = posterior$betaIn0 * posterior$Tin0 + posterior$betaIs0 * (posterior$Tinf - posterior$Tin0)
    ##posterior$Rt = posterior$betaInt * posterior$Tint + posterior$betaIst * (posterior$Tinf - posterior$Tint)

    posterior
}

fit.paramnames <- c("RIn0", "RInt", "R0", "Rt", "logHR", "HL", "DL",
                    "lockdownmort", "gamma.in0", "gamma.int")
keyparamnames <- c("betaIn0", "betaInt", "betaIs0", "betaIst", "R0", "Rt", "Tef0", "Teft",
                   "RIn0", "RIs0")
fitkeyparamnames <- c("R0", "Rt", "RIn0", "RInt", "gamma.in0", "gamma.int")

if (length(Es.time) > 0) {
    for (i in 1:length(Es.time)) {
        fit.paramnames <- c(fit.paramnames, paste("E", i, sep=""))
    }
}

simPlot <- function(labels, legendposition, th = th) {
    orig <- c(1.988216, 0.8245114, 3, 0.8, -3.152355, 7, 21.04874,
              3.377203, 0.3972412, 0.9936065)
    names(orig) <- fit.paramnames

    f <- function(state) { -diff(state$S) }
    
    range <- 600:1100
    df1 <- data.frame(x=range/10)
    df2 <- data.frame(x=range/10)
    df3 <- data.frame(x=range/10)
    df4 <- data.frame(x=range/10)
    df5 <- data.frame(x=range/10)
    df6 <- data.frame(x=range/10)
    
    ## 1
    ##  RIn0 = R0, RIs0 = 0
    ##  RIst = Rt, RIst = 0
    ##  Tin0 = Tint = 2.4
    p1 <- orig
    p1[1] = p1[3]
    p1[2] = p1[4]
    p1[9] = p1[10] = 1/2.4

    p1 <- transformParams(unlist(p1, use.names=FALSE))
    state <- calcNominalState(calculateModel(p1, 100))

    df1$i = f(state)[range]
    df1$Rt = state$Rt[range]
    df1$Re = state$Re[range]
    df1$hospi = state$hospi[range]
    df1$ci = N - state$S[range]
    
    ## 2
    ##  RIn0 = R0, RIs0 = 0
    ##  RIst = Rt, RIst = 0
    ##  betaIn0 = betaInt = 1
    p2 <- orig
    p2[1] = p2[3]
    p2[2] = p2[4]
    p2[9] = 1/2.4
    p2[10] = p2[1] / p2[2] * p2[9] # R0 * gamma0 = Rt * gammat

    p2 <- transformParams(unlist(p2, use.names=FALSE))
    state <- calcNominalState(calculateModel(p2, 100))

    i.jitter <- 10
    Rt.jitter <- 0.01
    
    df2$i = f(state)[range]
    df2$Rt = state$Rt[range]
    df2$Re = state$Re[range]
    df2$hospi = state$hospi[range]
    df2$ci = N - state$S[range]

    ## 3
    ##  RIn0 = R0, RIs0 = 0
    ##  RIst = Rt, RIst = 0
    ##  Tin0 = Tint = 2.4

    orig[4] = 0.6
    
    p3 <- orig
    p3[1] = p3[3]
    p3[2] = p3[4]
    p3[9] = p3[10] = 1/2.4

    p3 <- transformParams(unlist(p3, use.names=FALSE))
    state <- calcNominalState(calculateModel(p3, 100))

    df3$i = f(state)[range] + i.jitter
    df3$Rt = state$Rt[range] + Rt.jitter
    df3$Re = state$Re[range]
    df3$hospi = state$hospi[range]
    df3$ci = N - state$S[range]

    ## 4
    ##  RIn0 = R0, RIs0 = 0
    ##  RIst = Rt, RIst = 0
    ##  betaIn0 = betaInt = 1
    p4 <- orig
    p4[1] = p4[3]
    p4[2] = p4[4]
    p4[9] = 1/2.4
    p4[10] = p4[1] / p4[2] * p4[9] # R0 * gamma0 = Rt * gammat

    p4 <- transformParams(unlist(p4, use.names=FALSE))
    state <- calcNominalState(calculateModel(p4, 100))

    df4$i = f(state)[range] + i.jitter
    df4$Rt = state$Rt[range] + Rt.jitter
    df4$Re = state$Re[range]
    df4$hospi = state$hospi[range]
    df4$ci = N - state$S[range]

    ## 5
    ##  RIn0 = R0, RIs0 = 0
    ##  RIst = Rt, RIst = 0
    ##  Tin0 = Tint = 2.4

    orig[4] = 1.2
    
    p5 <- orig
    p5[1] = p5[3]
    p5[2] = p5[4]
    p5[9] = p5[10] = 1/2.4

    p5 <- transformParams(unlist(p5, use.names=FALSE))
    state <- calcNominalState(calculateModel(p5, 100))

    df5$i = f(state)[range] + i.jitter
    df5$Rt = state$Rt[range] + Rt.jitter
    df5$Re = state$Re[range]
    df5$hospi = state$hospi[range]
    df5$ci = N - state$S[range]

    ## 6
    ##  RIn0 = R0, RIs0 = 0
    ##  RIst = Rt, RIst = 0
    ##  betaIn0 = betaInt = 1
    p6 <- orig
    p6[1] = p6[3]
    p6[2] = p6[4]
    p6[9] = 1/2.4
    p6[10] = 1/7 ## p6[1] / p6[2] * p6[9] # R0 * gamma0 = Rt * gammat

    p6 <- transformParams(unlist(p6, use.names=FALSE))
    state <- calcNominalState(calculateModel(p6, 100))

    df6$i = f(state)[range] + i.jitter
    df6$Rt = state$Rt[range] + Rt.jitter
    df6$Re = state$Re[range]
    df6$hospi = state$hospi[range]
    df6$ci = N - state$S[range]

    top <- ggplot(df1,aes(x,i)) +
        geom_line(aes(linetype="i1", color="08")) +
        geom_line(data=df2, aes(linetype="i2", color="08")) +
        geom_line(data=df3, aes(linetype="i1", color="04")) +
        geom_line(data=df4, aes(linetype="i2", color="04")) +
        ##geom_line(data=df5, aes(linetype="i1", color="12")) +
        ##geom_line(data=df6, aes(linetype="i2", color="12")) +
        geom_vline(xintercept=(lockdown_offset + state$offset), linetype="dashed",
                   color='orange', size=0.5) +
        geom_vline(xintercept=(lockdown_offset + state$offset + lockdown_transition_period),
                   linetype="dashed", color='red', size=0.5) +
        scale_linetype_manual(values = c("i1" = "solid", "i2" = "dashed"), name='',
                              labels=expression(paste("Variable ",beta,", constant ",gamma),
                                                paste("Constant ",beta,", variable ",gamma))) +
        scale_color_manual(values = c("08" = "blue", "04" = "#33AA66", "12" = "#AA6633"), name='',
                           labels=c("Rt = 0.4", "Rt = 0.6")) +
        labs(linetype=NULL) +
        labs(color=expression(paste("Legend"))) +
        labs(y = "New infections per day") +
        theme(legend.position = legendposition) +
        theme(legend.title = element_blank())

    middle <- ggplot(df1,aes(x,hospi)) +
        geom_line(aes(linetype="i1", color="08")) +
        geom_line(data=df2, aes(linetype="i2", color="08")) +
        geom_line(data=df3, aes(linetype="i1", color="04")) +
        geom_line(data=df4, aes(linetype="i2", color="04")) +
        ##geom_line(data=df5, aes(linetype="i1", color="12")) +
        ##geom_line(data=df6, aes(linetype="i2", color="12")) +
        geom_vline(xintercept=(lockdown_offset + state$offset), linetype="dashed",
                   color='orange', size=0.5) +
        geom_vline(xintercept=(lockdown_offset + state$offset + lockdown_transition_period),
                   linetype="dashed", color='red', size=0.5) +
        theme(legend.position = "none") +
        labs(y = "Confirmed cases per day") +
        scale_linetype_manual(values = c("i1" = "solid", "i2" = "dashed"), name='',
                              labels=expression(paste("Variable ",beta,", constant ",gamma),
                                                paste("Constant ",beta,", variable ",gamma))) +
        scale_color_manual(values = c("08" = "blue", "04" = "#33AA66", "12" = "#AA6633"), name='',
                           labels=c("Rt = 0.4", "Rt = 0.8"))

    Rt <- ggplot(df1,aes(x,Rt)) +
        geom_line(aes(), color="blue") +
        geom_line(data=df3, aes(), color="#33AA66") +
        geom_vline(xintercept=(lockdown_offset + state$offset), linetype="dashed",
                   color='orange', size=0.5) +
        geom_vline(xintercept=(lockdown_offset + state$offset + lockdown_transition_period),
                   linetype="dashed", color='red', size=0.5) +
        scale_y_continuous(limits=c(0, 3.5)) +
        labs(y = "Rt") +
        labs(x = "Day")

    cowplot::plot_grid(top + theme(axis.text.x = element_blank(),
                                   axis.ticks.x = element_blank(),
                                   axis.title.x = element_blank()) + th,
                       middle + theme(axis.text.x = element_blank(),
                                      axis.ticks.x = element_blank(),
                                      axis.title.x = element_blank()) + th,
                       Rt + th,
                       nrow = 3, labels=labels, align="v", rel_heights=c(1,1,0.4))
 }

country2 = "BE"
source("../data/ecdc/data.R", chdir=T)

N <- N * 10

lockdown_offset <- 30
lockdown_transition_period <- 0

p1 <- simPlot("auto", "none", theme())

lockdown_offset <- 30 - 3.5
lockdown_transition_period <- 7

th <- theme(axis.title.y = element_blank())

p2 <- simPlot(c("d", "e", "f"), c(0.7, 0.7), th)

##pdf(paste("model-cmp.pdf", sep=""), width=8, height=8)
grid.arrange(p1, p2, nrow=1)
##dev.off()
