##
## Model for simulating exiting scenario's
## Assumes that model.R has been loaded already
##

## relax lockdown at this day
relax_offset <- as.numeric(as.Date("2020/6/30") - dstartdate) 
relax_E = 0.7 # 0 : no lockdown, 1 : lockdown 'lite' as now

## lift all measures at this day
lift_offset <- as.numeric(as.Date("2020/8/31") - dstartdate) 

calculateModel <- function(params, period)
{
    beta0 <- params[1]
    betat <- params[2]
    hosp_rate <- params[3]
    died_rate <- params[4]
    hosp_latency <- params[5]
    died_latency <- params[6]
    Tinf <- params[7]
    Tinc <- params[8]
    mort_lockdown_threshold <- params[9]
    hosp_rate_change <- params[10]

    a <- 1 / Tinc
    gamma <- 1 / Tinf

    hosp_cv_profile = calcConvolveProfile(-hosp_latency, 5)
    died_cv_profile = calcConvolveProfile(-died_latency, 5)

    padding = max(-hosp_cv_profile$kbegin, -died_cv_profile$kbegin) + 1

    state <- NULL
    state$S <- rep(N-Initial, padding + period)
    state$E <- rep(Initial, padding + period)
    state$I <- rep(0, padding + period)
    state$R <- rep(0, padding + period)
    state$hosp <- rep(0, padding + period)
    state$died <- rep(0, padding + period)
    state$i <- padding + 1

    data_offset = InvalidDataOffset

    beta = beta0	
    for (i in (padding + 1):(padding + period)) {
        betai = i - data_offset - lockdown_offset
	if (betai < 0) {
            beta = beta0
	} else {
            if (betai >= lockdown_transition_period) {
                if (i - data_offset > relax_offset) {
                    if (i - data_offset > lift_offset) {
                        beta = beta0
                    } else {
                        beta = max(0.001, relax_E * betat + (1 - relax_E) * beta0)
                    }
                } else {
                    beta = betat
                }
            } else {
                fract0 = (lockdown_transition_period - betai) / lockdown_transition_period
                fractt = betai / lockdown_transition_period
                beta = fract0 * beta0 + fractt * betat
            }
	}

	state <- simstep(state, 1, N, beta, a, gamma)

    	s = convolute(state$S, i, hosp_cv_profile)
	state$hosp[i] <- (N - s) * hosp_rate

	r = convolute(state$R, i, died_cv_profile)
	state$died[i] <- r * died_rate

## assuming lockdown decision was based on a cumulative mort count, with some
## uncertainty on the exact value due to observed cumulative mort count being a
## sample
	if (data_offset == InvalidDataOffset && state$died[i] > mort_lockdown_threshold) {
            data_offset = i - lockdown_offset
	}
    }

    state$deadi <- numeric(length(state$died))
    state$hospi <- numeric(length(state$hosp))

    for (i in (padding + 1):(padding + period)) {
    	state$deadi[i] <- if (i == 1) state$died[1] else state$died[i] - state$died[i - 1]
	state$hospi[i] <- if (i == 1) state$hosp[1] else state$hosp[i] - state$hosp[i - 1]
    }

    cumDiff <- 0
    inv_hosp_rate_change <- 1 / hosp_rate_change
    for (i in (round(data_offset + 12)):length(state$hospi)) {
    	nhospi <- state$hospi[i] * inv_hosp_rate_change
	diff <- nhospi - state$hospi[i]
    	state$hospi[i] <- nhospi
	cumDiff <- cumDiff + diff
	state$hosp[i] <- state$hosp[i] + cumDiff
    }

    state$padding <- padding
    state$offset <- data_offset

    state
}
