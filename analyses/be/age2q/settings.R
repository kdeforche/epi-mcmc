## R dir (all scripts)
Rdir <- "../../../R/"

##
HospLabel <- "Cases"

## Model to be used for data fitting (fitMCMC.R, evalMCMC.R)
fitmodel <- "../../../R/models/model-age-groups2q.R"

## Data to be used (fitMCMC.R, evalMCMC.R)
data <- "../../../R/data/be/data-age-all.R"

## Output file for MCMC samples (fitMCMC.R)
outputfile <- paste("BE_run.csv", sep='')

## Input file for MCMC samples (evalMCMC.R, predictMCMC.R)
#inputfiles <- c("run4.csv", "1/run4.csv", "2/run4.csv", "3/run4.csv", "4/run4.csv", "5/run4.csv", "6/run4.csv")
inputfiles <- c(outputfile)

## Output directory for graphs, etc... (used by eval.R and predict.R)
outputdir <- paste("BE_results", sep='')

## Sample size used for creating density plots (evalMCMC.R)
densityPlotSampleSize <- 500

## Sample size used for creating density plots (predictMCMC.R)
quantilePlotSampleSize <- 1000

ifr_epsilon = 0.12

FitTotalPeriod = 300
case_nbinom_size1 = 0.01
case_nbinom_sizey2 = 1
case_nbinom_sizeo2 = 20
hosp_nbinom_size1 = 20 
hosp_nbinom_size2 = 40 
mort_nbinom_size = 120

if (.Device == "null device") {
    x11(width=15, height=8)
}
