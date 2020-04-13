source("data.R")
source("model.R")

options(scipen=999)

calcloglML <- function(params) {
    hosp_rate <- params[3]
    died_rate <- params[4]
    hosp_latency <- params[5]
    died_latency <- params[6]

    Tinf <- 7
    Tinc <- 5
    mort_lockdown_threshold = 5
    hosp_rate_change = 1.1

    beta0 <- params[1] / Tinf
    betat <- params[2] / Tinf

    logl <- calclogl(c(beta0, betat, hosp_rate, died_rate, hosp_latency, died_latency, Tinf, Tinc, mort_lockdown_threshold, hosp_rate_change))

    logl
}

o <- NULL
o$par <- c(4, 0.7, 0.03, 0.01, 10, 11)
calcloglML(o$par)

control <- NULL
control$maxit <- 10000
control$fnscale <- -1
control$parscale <- c(0.1, 0.1, 1/100, 1/300, 1, 1)
o <- optim(o$par, calcloglML, control=control)
o
calcloglML(o$par)
graphs()
