calcbeta.orig <- calcbeta

relax.start <- as.numeric(as.Date("2020/5/1") - dstartdate)
relax.measures.E <- 0.80 # reduce to 80%
relax.end <- as.numeric(as.Date("2020/8/1") - dstartdate)

calcbeta.relax <- function(time, beta0, betat)
{
    if (time < relax.start)
        return(calcbeta.orig(time, beta0, betat))

    if (time > relax.end)
        return(beta0)
    
    beta = max(0.001, relax.measures.E * betat + (1 - relax.measures.E) * beta0)

    beta
}
