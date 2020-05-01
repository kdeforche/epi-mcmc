# epi-mcmc
Scripts for estimating and visualizing epidemiological models of an epidemic (developed for COVID-19)

## Background

These scripts were developed to get insight in the evolution of the
Belgian COVID-19 epidemic during March-April 2020. The description of
the method and results for the Belgian COVID-19 epidemic are described
here: https://www.genomedetective.com/covid19-epi/be/report.pdf

Since the epidemic is still ongoing, the models may be expanded with
new parameters to model social distancing changes as restrictions are
being lifted (and/or resset). Visualizations are still evolving too.

## Description

The models are estimated using MCMC from incidence data of
hospitalization or case count (as a time-shifted proxy for infections)
and death counts. They have been designed to model the introduction of
social distancing ("lockdown" or other) and consider a different R0
before and after the introduction of the measures.

The model uses 'hospitalization' as name for the quantity that is
either hospitalization or positive test.

Three different model flavours are currently available.

### Simple model

This model is implemented as 'R/models/model.R'

The following are the parameters estimated by the model. The names
reflect commonly used terminology in SEIR models.

| Variable | Description                                       |
|:--------:|---------------------------------------------------|
| beta0    | beta (daily reproduction number) before lockdown  |
| betat    | beta (daily reproduction number) after lockdown   |
| logHR    | log(hospitalization rate)                         |
| logHRDR  | log(IFR / hospitalization rate)                   |
| HL       | hospitalization latency                           |
| DL       | death latency                                     |
| Tinf     | infectious period                                 |
| Tinc     | latent period                                     |
| lockdownmort | cumulative deaths to start lockdown           |

The following are derived parameters (some priors are based on these,
actually)

| Variable | Description                                       |
|:--------:|---------------------------------------------------|
| R0       | basic reproduction number before lockdown         |
| Rt       | basic reproduction number after lockdown          |
| HR       | hospitalization rate                              |
| IFR      | infection fatality rate (= 'DR')                  |
| Et       | Rt / R0                                           |

The data to fit the simple model needed are shown for example in 'R/data/de/data.R':

* *dhospi* : vector of incidence of hospitalizations
* *dhosp* : cumsum(dhospi)
* *dmorti* : vector of incidence of deaths
* *dmort* : cumsum(dmorti)
* *dstartdate* : first day of the above data vectors
* *N* : population size
* *lockdown_offset* : index (since dstartdate) of start of lockdown transition period
* *lockdown_transition_period* : duration of lockdown transition period
* *total_deaths_at_lockdown* : deaths at lockdown (dmort[lockdown_offset]), under the assumption that deaths were the reason for the government to take measures (WHO recommendation)

If you wish to setup an analysis for a new data set, then find a good
representative data set. Sometimes the national health agency provides
good data, our otherwise it is our experience that Wikipedia is a good
source of the data (go to SOURCE to scrape the data that is shown in
the graphs). For setting lockdown_offset and
lockdown_transition_period, the Google Mobility Reports are a fairly
reliable data source.

### Simple model with effective population size

This model is implemented as 'R/models/model-Ne.R'

This model adds a single parameter to co-estimate an effective
population size. This is typically needed when fitting the model to
case data which reaches a peak because of immunity, since an SEIR
model assumes a homogeneous mixing of the population. In practice,
this is never the case.

The effective population size can be lower than the real population
size and changes the rate at which infectious people are infecting
susceptible people by considering a lower 'relevant' total population
for both groups.

| Variable | Description                                       |
|:--------:|---------------------------------------------------|
| Nef      | Factor (0 - 1) of Ne relative to N                |

An example output of this model, fitting the Swedish data, is briefly described here:
https://twitter.com/houterkabouter/status/1256123165346062336?s=20

### Age-structured model

This model is implemented as 'R/models/model-age.R'

This is the model that has been described in the paper linked
above. It considers two age groups separately: 'y' (younger) versus
'o' (older).

| Variable | Description                                                               |
|:--------:|---------------------------------------------------------------------------|
| y.beta0  | beta (daily reproduction number) before lockdown, younger group           |
| y.betat  | beta (daily reproduction number) after lockdown, younger group            |
| o.beta0  | beta (daily reproduction number) before lockdown, older group             |
| o.betat  | beta (daily reproduction number) after lockdown, older group              |
| yo.beta0 | beta (daily reproduction number) before lockdown, younger to older group  |
| yo.betat | beta (daily reproduction number) after lockdown, younger to older group   |
| logHL    | log(hospitalization rate), younger group                                  |
| logHRyDR | log(IFR / younger hospitalization rate), younger group                    |
| logHRoDR | log(IFR / younger hospitalization rate), older group                      |
| y.HL     | hospitalization latency, younger group                                    |
| y.DL     | death latency, younger group                                              |
| o.HL     | hospitalization latency, older group                                      |
| o.DL     | death latency, older group                                                |
| Tinf     | infectious period                                                         |
| Tinc     | latent period                                                             |
| lockdownmort | cumulative deaths to start lockdown                                   |
| logfHo   | log(older hospitalization rate / younger hospitalization rate)            |

It requires separate data sets of death and hospitalization incidence
for each group (y.*x* or o.*x*) :

* *y.dhospi*, *o.dhospi*
* *y.dhosp*, *o.dhosp*
* *y.dmorti*, *o.dmorti*
* *y.dmort*, *o.dmort*
* *y.N*, *o.N*

## Usage

### Estimating a model

To estimate the model from data, you will need to follow these steps:

* Setup a new dataset file (add to 'R/data/' ) or use an existing one from that folder
* Create and edit a 'settings.R' file (copy from an existing one in the analysis folder)
* Rscript R/MCMC/fitMCMC.R

This will start to sample from the posterior distribution. You will
need to collect a sufficient number before you can start making useful
conclusions. This takes multiple hours (and running a few in parallel
on a server is a good idea). If you are new to MCMC, read a bit about
MCMC and diagnostics while waiting ;-)

### Visualizations

Visualizations are less well organized (and seemingly in constant flux). There are two
main scripts (but they will not run entirely out of the box and have a lot of things in them
that can be commented in and out).

* 'R/MCMC/evalMCMC.R' : basic MCMC diagnostics. Use this first to monitor the
  progress of MCMC and visualize posterior estimate histograms.

* 'R/MCMC/predictMCMC.R' : creates the plots that are also in the publication

