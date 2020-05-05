## R dir (all scripts)
Rdir <- "../../../R/"

##
HospLabel <- "New cases"

## Model to be used for data fitting (fitMCMC.R, evalMCMC.R)
fitmodel <- "../../../R/models/model-Ne2-inf.R"

## Data to be used (fitMCMC.R, evalMCMC.R)
data <- "../../../R/data/nl/data.R"

## Output file for MCMC samples (fitMCMC.R)
outputfile <- "ne2-inf.csv"

## Input file for MCMC samples (evalMCMC.R, predictMCMC.R)
#inputfiles <- c("run4.csv", "1/run4.csv", "2/run4.csv", "3/run4.csv", "4/run4.csv", "5/run4.csv", "6/run4.csv")
inputfiles <- c(outputfile)

## Output directory for graphs, etc... (used by eval.R and predict.R)
outputdir <- "results"

## Sample size used for creating density plots (evalMCMC.R)
densityPlotSampleSize <- 500

## Sample size used for creating density plots (predictMCMC.R)
quantilePlotSampleSize <- 50

FitTotalPeriod <- 130

hosp_nbinom_size = 15
mort_nbinom_size = 30

if (.Device == "null device") {
    x11(width=15, height=12)
}
