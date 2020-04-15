calcbetas.age.orig <- calcbetas.age

relax.start <- as.numeric(as.Date("2020/5/1") - dstartdate)
relax.measures.y.E <- 0.5 # reduce to 50%
relax.measures.o.E <- 1   # keep 100%
relax.measures.yo.E <- 1  # keep 100%
relax.end <- as.numeric(as.Date("2020/8/1") - dstartdate)

calcbetas.age.relax <- function(time, betay0, betayt, betao0, betaot, betayo0, betayot)
{
    if (time < relax.start)
        return(calcbetas.age.orig(time, betay0, betayt, betao0, betaot, betayo0, betayot))

    if (time > relax.end)
        c(betay0, betao0, betayo0)
    
    y.beta = max(0.001, relax.measures.y.E * betayt + (1 - relax.measures.y.E) * betay0)
    o.beta = max(0.001, relax.measures.o.E * betaot + (1 - relax.measures.o.E) * betao0)
    yo.beta = max(0.001, relax.measures.yo.E * betayot + (1 - relax.measures.yo.E) * betayo0)

    c(y.beta, o.beta, yo.beta)
}
