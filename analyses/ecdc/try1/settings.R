## Country
args <- commandArgs(TRUE)
country2 <- args[length(args)]
print(country2)

## R dir (all scripts)
Rdir <- "../../../R/"

##
HospLabel <- "New cases"

## Model to be used for data fitting (fitMCMC.R, evalMCMC.R)
fitmodel <- "../../../R/models/model-Ne3-inf.R"

## Data to be used (fitMCMC.R, evalMCMC.R)
data <- "../../../R/data/ecdc/data.R"

## Output file for MCMC samples (fitMCMC.R)
outputfile <- paste(country2, "_run.csv", sep='')

## Input file for MCMC samples (evalMCMC.R, predictMCMC.R)
#inputfiles <- c("run4.csv", "1/run4.csv", "2/run4.csv", "3/run4.csv", "4/run4.csv", "5/run4.csv", "6/run4.csv")
inputfiles <- c(outputfile)

## Output directory for graphs, etc... (used by eval.R and predict.R)
outputdir <- paste(country2, "_results", sep='')

## Sample size used for creating density plots (evalMCMC.R)
densityPlotSampleSize <- 500

## Sample size used for creating density plots (predictMCMC.R)
quantilePlotSampleSize <- 50

FitTotalPeriod <- 160

hosp_nbinom_size = 5 
mort_nbinom_size = 25

if (.Device == "null device") {
    x11(width=15, height=12)
}
