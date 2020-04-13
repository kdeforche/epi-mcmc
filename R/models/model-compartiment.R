N <- 11.5E6

j.N <- 8.55E6
o.N <- N - j.N

# factor for died_rate with respect to estimated died_rate
j.died_rate_factor = 0.09412495/0.66
o.died_rate_factor = 2.34371716/0.66

o.beta0.factor <- 0.7
j.beta0.factor <- (N - o.beta0.factor * o.N) / j.N

o.betat.factor <- 1
j.betat.factor <- (N - o.betat.factor * o.N) / j.N

# factor for hosp_rate with respect to estimated hosp_rate
# assuming here also that same fraction of hosptilazations for deaths, but maybe
# younger people have higher hospitalization / dead ratio than elder ?
j.hosp_rate_factor = 0.09/0.66
o.hosp_rate_factor = 2.34/0.66

# measure efficiency : ratio of efficiency compare to lockdown efficiency
# lockdown efficiency = 1 - betat / beta0
# 1: betat; 0: beta0
o.E <- 1.0
j.E <- 0.1
jo.E <- 1.0

# let's take 13/03 as lockdown date (we're doing it gradual anyway)
lockdown_offset <- 3
InvalidDataOffset <- 10000
Initial <- 1

# date at which we start with this
start_offset <- as.numeric(as.Date("2020/4/30") - as.Date("2020/3/10"))

# date at which we go out of lockdown
end_offset <- as.numeric(as.Date("2020/9/30") - as.Date("2020/3/10"))

calcConvolveProfile <- function(latency, latency_sd)
{
    kbegin = floor(latency - latency_sd * 1.5)
    kend = min(0, ceiling(latency + latency_sd * 1.5))

    result = NULL
    result$kbegin = kbegin
    result$kend = kend
    result$values = numeric(kend - kbegin)
    i = 1
    for (k in kbegin:kend) {
       result$values[i] = pnorm(k-0.5, mean=latency, sd=latency_sd) - pnorm(k+0.5, mean=latency, sd=latency_sd)
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

simstep <- function(state, j.beta, o.beta, jo.beta, a, gamma)
{
  i = state$i

  calcjN = state$j.S[i] + state$j.I[i] + state$j.E[i] + state$j.R[i]
  calcoN = state$o.S[i] + state$o.I[i] + state$o.E[i] + state$o.R[i]

  j.got_infected = j.beta * state$j.I[i] * state$j.S[i] / calcjN
  j.got_infectious = a * state$j.E[i]
  j.got_removed = gamma * state$j.I[i]

  j.deltaS = -j.got_infected
  j.deltaE = j.got_infected - j.got_infectious
  j.deltaI = j.got_infectious - j.got_removed
  j.deltaR = j.got_removed

  state$j.S[i + 1] = state$j.S[i] + j.deltaS
  state$j.E[i + 1] = state$j.E[i] + j.deltaE
  state$j.I[i + 1] = state$j.I[i] + j.deltaI
  state$j.R[i + 1] = state$j.R[i] + j.deltaR

  if (calcoN > 0) {
    o.got_infected = (o.beta * state$o.I[i] + jo.beta * state$j.I[i]) * state$o.S[i] / calcoN
    o.got_infectious = a * state$o.E[i]
    o.got_removed = gamma * state$o.I[i]

    o.deltaS = -o.got_infected
    o.deltaE = o.got_infected - o.got_infectious
    o.deltaI = o.got_infectious - o.got_removed
    o.deltaR = o.got_removed

    state$o.S[i + 1] = state$o.S[i] + o.deltaS
    state$o.E[i + 1] = state$o.E[i] + o.deltaE
    state$o.I[i + 1] = state$o.I[i] + o.deltaI
    state$o.R[i + 1] = state$o.R[i] + o.deltaR
  }

  i = i + 1
  state$i = i

  state
}

calculateModel <- function(E0, params, period)
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

    padding = floor(max(hosp_latency, died_latency) * 3)

    state <- NULL

    state$j.S <- rep(j.N - E0, padding + period)
    state$j.E <- rep(E0/2, padding + period)
    state$j.I <- rep(0, padding + period)
    state$j.R <- rep(0, padding + period)
    state$j.hosp <- rep(0, padding + period)
    state$j.died <- rep(0, padding + period)

    state$o.S <- rep(o.N - E0/2, padding + period)
    state$o.E <- rep(E0/2, padding + period)
    state$o.I <- rep(0, padding + period)
    state$o.R <- rep(0, padding + period)
    state$o.hosp <- rep(0, padding + period)
    state$o.died <- rep(0, padding + period)

    state$i <- padding + 1

    beta_transition = 7
    data_offset = InvalidDataOffset

    j.beta = beta0
    o.beta = beta0 * o.beta0.factor
    jo.beta = 0

    j.hr = j.hosp_rate_factor * hosp_rate
    j.dr = j.died_rate_factor * died_rate

    o.hr = o.hosp_rate_factor * hosp_rate
    o.dr = o.died_rate_factor * died_rate

    for (i in (padding + 1):(padding + period)) {
        betai = i - data_offset - lockdown_offset

        if (betai < 0) {
	   j.beta = beta0 * j.beta0.factor
	   o.beta = beta0 * o.beta0.factor
	} else {
	   if (betai >= beta_transition) {
	      j.beta = betat * j.betat.factor
	      o.beta = betat * o.betat.factor
   	   } else {
	      fract0 = (beta_transition - betai) / beta_transition
	      fractt = betai / beta_transition
	      j.beta = fract0 * beta0 * j.beta0.factor + fractt * betat * j.betat.factor
	      o.beta = fract0 * beta0 * o.beta0.factor + fractt * betat * o.betat.factor
	   }
	}

	if (i - data_offset >= start_offset) {
	   if (i - data_offset < end_offset) {
	      j.beta = (1 - j.E) * beta0 * j.beta0.factor + j.E * betat * j.betat.factor
	      o.beta = (1 - o.E) * beta0 * o.beta0.factor + o.E * betat * o.betat.factor
	      jo.beta = (1 - jo.E) * beta0 * o.beta0.factor + jo.E * betat * o.betat.factor
          } else {
             # we should remove the compartments, for now let's assume that
	     # this is a good approximation:
	     j.beta = beta0 * j.beta0.factor
	     o.beta = beta0 * o.beta0.factor
	     jo.beta = beta0 * o.beta0.factor
          }
	}

	state <- simstep(state, j.beta, o.beta, jo.beta, a, gamma)

    	s = convolute(state$j.S, i, hosp_cv_profile)
	state$j.hosp[i] <- (j.N - s) * j.hr

	r = convolute(state$j.R, i, died_cv_profile)
	state$j.died[i] <- r * j.dr

	# if (state$j.died[i] < state$j.died[i - 1]) {
	#    print(c("===================== j.died[i] < j.died[i - 1]"))
	#    print(params)
	# }

    	s = convolute(state$o.S, i, hosp_cv_profile)
	state$o.hosp[i] <- (o.N - s) * o.hr

	r = convolute(state$o.R, i, died_cv_profile)
	state$o.died[i] <- r * o.dr

	# if (state$o.died[i] < state$o.died[i - 1]) {
	#    print(c("===================== o.died[i] < o.died[i - 1]"))
	#    print(state$o.R)
	#    print(state$o.died)
	#    print(params)
	# }

	# assuming lockdown decision was based on a cumulative mort count, with some
	# uncertainty on the exact value due to observed cumulative mort count being a
	# sample
	if (data_offset == InvalidDataOffset && state$j.died[i] > mort_lockdown_threshold) {
	   data_offset = i - lockdown_offset
	}
    }

    state$j.deadi <- numeric(length(state$j.died))
    state$j.hospi <- numeric(length(state$j.hosp))

    for (i in (padding + 1):(padding + period)) {
    	state$j.deadi[i] <- if (i == 1) state$j.died[1] else
			 state$j.died[i] - state$j.died[i - 1]
	state$j.hospi[i] <- if (i == 1) state$j.hosp[1] else
			 state$j.hosp[i] - state$j.hosp[i - 1]
    }

    state$o.deadi <- numeric(length(state$o.died))
    state$o.hospi <- numeric(length(state$o.hosp))

    for (i in (padding + 1):(padding + period)) {
    	state$o.deadi[i] <- if (i == 1) state$o.died[1] else
			 state$o.died[i] - state$o.died[i - 1]
	state$o.hospi[i] <- if (i == 1) state$o.hosp[1] else
			 state$o.hosp[i] - state$o.hosp[i - 1]
    }

    cumDiff <- 0
    inv_hosp_rate_change <- 1 / hosp_rate_change
    for (i in (round(data_offset + 12)):length(state$j.hospi)) {
    	nhospi <- state$j.hospi[i] * inv_hosp_rate_change
	diff <- nhospi - state$j.hospi[i]
    	state$j.hospi[i] <- nhospi
	cumDiff <- cumDiff + diff
	state$j.hosp[i] <- state$j.hosp[i] + cumDiff
    }

    for (i in (round(data_offset + 12)):length(state$o.hospi)) {
    	nhospi <- state$o.hospi[i] * inv_hosp_rate_change
	diff <- nhospi - state$o.hospi[i]
    	state$o.hospi[i] <- nhospi
	cumDiff <- cumDiff + diff
	state$o.hosp[i] <- state$o.hosp[i] + cumDiff
    }

    state$padding <- padding
    state$offset <- data_offset

    state
}

transformParams <- function(params)
{
    result = params
    result[4] = exp(params[3] + params[4])
    result[3] = exp(params[3])
    result[10] = exp(params[10])

    result
}
