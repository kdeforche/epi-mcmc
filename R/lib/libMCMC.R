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
    print(summary(posterior))

    posterior <- invTransformParams(posterior)
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
        state <- calcNominalState(calculateModel(params, 200))
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
        state <- calcNominalState(calculateModel(params, 200))
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

grey95 <- "grey70"
grey75 <- "grey60"

makePlot <- function(sample, dateRange, fun, colour, titles, date_markers, lty)
{
    period <- as.numeric(dateRange[2] - dateRange[1])
    qd <- data.frame(quantileData(sample, fun, period, c(0.05, 0.25, 0.5, 0.75, 0.95)))

    colnames(qd) <- c("q5", "q25", "q50", "q75", "q95")
    qd$x <- seq(dateRange[1], dateRange[1] + dim(qd)[1] - 1, 1)

    result <- ggplot(qd) + aes(x = x) +
        geom_ribbon(aes(ymin = q5, ymax=q95, fill=grey95), alpha=0.4) +
        geom_ribbon(aes(ymin = q25, ymax=q75, fill=grey75), alpha=0.4) +
        geom_line(aes(y = q50), linetype=lty, size = 1, colour=colour) +
        labs(x="Date") +
        labs(y=titles[1]) +
        ggtitle(titles[2]) +
        scale_fill_identity(name='Uncertainty', guide=FALSE, labels=c('90%', '50%')) +
        scale_x_date(date_breaks="months", date_labels = "%e %b") +
        theme(legend.position = c(0.8, 0.85)) +
        guides(linetype=guide_legend(keywidth = 3, keyheight = 1))

    for (i in 1:length(date_markers$pos)) {
        pos = date_markers$pos[i]
        color = date_markers$color[i]
        result = result + geom_vline(xintercept=pos, linetype="dashed",
                                     color=color, size=1.3)
    }

    result
}

addExtraPlot <- function(plot, sample, dateRange, fun, colour, lty)
{
    period <- as.numeric(dateRange[2] - dateRange[1])
    qd <- data.frame(quantileData(sample, fun, period, c(0.5)))

    colnames(qd) <- c("q50")
    
    result <- plot +
        geom_line(aes(y = qd$q50, linetype = lty), colour=colour, size = 1)

    result
}

addExtraPlotQ <- function(plot, sample, dateRange, fun, colour, lty)
{
    period <- as.numeric(dateRange[2] - dateRange[1])
    qd <- data.frame(quantileData(sample, fun, period, c(0.05, 0.25, 0.5, 0.75, 0.95)))

    colnames(qd) <- c("q5", "q25", "q50", "q75", "q95")
    
    result <- plot +
        geom_ribbon(aes(ymin = qd$q5, ymax=qd$q95, fill=grey95), alpha=0.4) +
        geom_ribbon(aes(ymin = qd$q25, ymax=qd$q75, fill=grey75), alpha=0.4) +
        geom_line(aes(y = qd$q50), linetype = lty, colour=colour, size = 1)

    result
}

addExtraPlotQ2 <- function(plot, sample, dateRange, fun, colour, lty)
{
    period <- as.numeric(dateRange[2] - dateRange[1])
    qd <- data.frame(quantileData(sample, fun, period, c(0.05, 0.25, 0.5, 0.75, 0.95)))

    colnames(qd) <- c("cq5", "cq25", "cq50", "cq75", "cq95")
    
    result <- plot +
        geom_ribbon(aes(ymin = qd$cq5, ymax=qd$cq95, fill="grey70"), alpha=0.4) +
        geom_ribbon(aes(ymin = qd$cq25, ymax=qd$cq75, fill="grey60"), alpha=0.4) +
        geom_line(aes(y = qd$cq50), linetype = lty, colour=colour, size = 1)

    result
}
