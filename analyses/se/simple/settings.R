## R dir (all scripts)
Rdir <- "../../../../R/"

##
HospLabel <- "New cases"

## Model to be used for data fitting (fitMCMC.R, evalMCMC.R)
fitmodel <- "../../../../R/models/model.R"

## Data to be used (fitMCMC.R, evalMCMC.R)
data <- "../../../../R/data/se/data.R"

## Output file for MCMC samples (fitMCMC.R)
outputfile <- "run3.csv"

## Input file for MCMC samples (evalMCMC.R, predictMCMC.R)
inputfiles <- c("run3.csv")

## Output directory for graphs, etc... (used by eval.R and predict.R)
outputdir <- "results"

## Sample size used for creating density plots (evalMCMC.R)
densityPlotSampleSize <- 500

## Sample size used for creating density plots (predictMCMC.R)
quantilePlotSampleSize <- 50

x11(width=15, height=12)
