##
## Reads and combines one or MCMC sample files
## 
readData <- function(files)
{
    posterior <- read.csv(files[1])
    print(dim(posterior))
    c = dim(posterior)[1]

    posterior <- posterior[(c/8):c,]

    if (length(files) > 1) {
        for (i in 2:length(files)) {
            posteriorb <- read.csv(files[i])
            print(dim(posteriorb))
            c = dim(posteriorb)[1]

            posteriorb <- posteriorb[(c/8):c,]
            posterior <- rbind(posterior, posteriorb)
        }
    }
    summary(posterior)

    posterior$y.R0 = posterior$betay0 * posterior$Tinf
    posterior$y.Rt = posterior$betayt * posterior$Tinf
    posterior$o.R0 = posterior$betao0 * posterior$Tinf
    posterior$o.Rt = posterior$betaot * posterior$Tinf
    posterior$yo.R0 = posterior$betayo0 * posterior$Tinf
    posterior$yo.Rt = posterior$betayot * posterior$Tinf
    posterior$IFR = exp(posterior$logHR + posterior$logHRDR)
    posterior$HR = exp(posterior$logHR)
    posterior$y.Et = posterior$y.Rt / posterior$y.R0
    posterior$o.Et = posterior$o.Rt / posterior$o.R0
    posterior$yo.Et = posterior$yo.Rt / posterior$yo.R0
    posterior
}

takeAndPad <- function(data, offset, l)
{
    result <- data[offset:(min(l,length(data)))]
    if (length(result != l))
        result <- c(result, rep(NA, l - length(result)))
    result
}

##
## Draws a single model curve for incidence of hospitalisations and/or deaths
## 
graph_daily <- function(state, first, color, last) {
    days <- seq(dstartdate, as.Date("2020/5/1"), 1)

    period <- length(state$hospi)
    len <- period - state$offset + 1

    hosp_color = color
    mort_color = NULL
    if (is.null(color)) {
        hosp_color = alpha("red", 0.1)
        mort_color = alpha("blue", 0.1)
    }

    if (first) {
        plot(days[1:len],state$hospi[state$offset:period], type='l', col=hosp_color,
	     xlab='Date', ylab='Count', ylim=c(0, 800))
    } else {
        lines(days[1:len],state$hospi[state$offset:period], type='l', col=hosp_color)
        ## lines(days[1:len],(state$E[state$offset:period] + state$I[state$offset:period])/N,
        ##       type='l', col=alpha("gray", 0.5))
    }

    if (!is.null(mort_color)) {
        lines(days[1:len],state$deadi[state$offset:period], type='l', col=mort_color)

    	if (last) {
            lines(days[1:length(dmorti)],dmorti,col=c("black"),type='o')
        }
    }
    if (last) {
        lines(days[1:length(dhospi)],dhospi,col=c("black"),type='o')
    }
}

##
## Draws a single model curve for cumulative deaths
##
graph_cum <- function(state, first, color, last) {
    days <- seq(dstartdate, as.Date("2020/8/1"), 1)

    period <- length(state$hospi)
    len <- period - state$offset + 1

    mort_color = NULL
    if (is.null(color)) {
        mort_color = alpha("blue", 0.1)
    }

    if (first) {
        plot(days[1:len],state$died[state$offset:period], type='l', col=mort_color,
	     xlab='Date', ylab='Count', ylim=c(0, 100000))
    } else {
        lines(days[1:len],state$died[state$offset:period], type='l', col=mort_color)
    }

    if (last) {
        lines(days[1:length(dmort)],dmort,col=c("black"),type='o')
    }
}

##
## Draws daily incidence death for a posterior sample of parameter estimates
##
predict_daily_plot <- function(sample, overlay, color) {
    first <- !overlay

    for (i in 1:(dim(sample)[1])) {
        params <- sample[i,]
        params <- transformParams(unlist(params, use.names=FALSE))
        state <- calculateModel(params, 200)
        state$hosp <- state$y.hosp + state$o.hosp
        state$died <- state$y.died + state$o.died
        state$hospi <- state$y.hospi + state$o.hospi
        state$deadi <- state$y.deadi + state$o.deadi
        graph_daily(state, first, color, i == dim(sample)[1])
        first <- FALSE
    }

    if (!overlay) {
        title(paste(c(HospLabel, 'and deaths per day')))
        legend("topleft", inset=0.02, legend=c(HospLabel, "Deaths"),
               col=c("red", "blue"),lty=1)
    }
}

##
## Draws cumulative death for a posterior sample of parameter estimates
##
predict_cum_plot <- function(sample, overlay, color) {
    first <- !overlay

    ests <- c()
    for (i in 1:(dim(sample)[1])) {
        params <- sample[i,]
        params <- transformParams(unlist(params, use.names=FALSE))
        state <- calculateModel(params, 200)
        ests <- c(ests, state$died[length(state$died - 1)])
        graph_cum(state, first, color, i == dim(sample)[1])
        first <- FALSE
    }
    ests
}

##
## Calculates quantiles of the given fun, over the given period, for the posterior state sample
##
quantileData <- function(sample, fun, period, quantiles)
{
    simperiod = period * 3
    
    sampleSize <- dim(sample)[1]
    data <- matrix(nrow=simperiod, ncol=sampleSize)

    maxOffset <- 0
    for (i in 1:sampleSize) {
        params <- sample[i,]
        params <- transformParams(unlist(params, use.names=FALSE))
        state <- calculateModel(params, simperiod)
        v <- takeAndPad(fun(state), state$offset, simperiod)
        maxOffset <- max(maxOffset, state$offset)
        data[,i] = v
    }

    result <- data.frame()
    result <- matrix(nrow=period, ncol=length(quantiles))
    for (j in 1:(dim(result)[1])) {
        result[j,] = quantile(data[j,], quantiles, na.rm=T)
    }
    result
}

##
## Creates a banded plot of the given fun, over the given period, for te posterior state sample 
##
makePlot <- function(sample, dateRange, fun, colour, titles, date_markers)
{
    period <- as.numeric(dateRange[2] - dateRange[1])
    qd <- data.frame(quantileData(sample, fun, period, c(0.05, 0.25, 0.5, 0.75, 0.95)))

    colnames(qd) <- c("q5", "q25", "q50", "q75", "q95")
    qd$x <- seq(dateRange[1], dateRange[1] + dim(qd)[1] - 1, 1)
    
    result <- ggplot(qd) + aes(x = x) +
        geom_ribbon(aes(ymin = q5, ymax=q95), alpha=0.4, fill="grey70") +
        geom_ribbon(aes(ymin = q25, ymax=q75), alpha=0.4, fill="grey60") +
        geom_line(aes(y = q50), colour=colour, size = 1) +
        labs(x="Date") +
        labs(y=titles[1]) +
        ggtitle(titles[2]) +
        scale_x_date(date_breaks="months", date_labels = "%e %b")

    for (i in 1:length(date_markers$pos)) {
        pos = date_markers$pos[i]
        color = date_markers$color[i]
        result = result + geom_vline(xintercept=pos, linetype="dashed",
                                     color=color, size=1.3)
    }

    result
}
