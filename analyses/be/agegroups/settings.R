## R dir (all scripts)
Rdir <- "../../../R/"

##
HospLabel <- "New hospitalisations"

## Model to be used for data fitting (fitMCMC.R, evalMCMC.R)
fitmodel <- paste(Rdir, "models/model-age-groups.R", sep="")

## Data to be used (fitMCMC.R, evalMCMC.R)
data <- paste(Rdir, "data/be/data-age.R", sep="")

## Output file for MCMC samples (fitMCMC.R)
outputfile <- "run.csv"
truncate <- T # whether the current file should be appended or instead truncated

## Input file for MCMC samples (evalMCMC.R, predictMCMC.R)
#inputfiles <- c(outputfile, "1/run.csv", "2/run.csv", "3/run.csv", "4/run.csv", "5/run.csv", "6/run.csv")
inputfiles <- c(outputfile)

## Output directory for graphs, etc... (used by eval.R and predict.R)
outputdir <- "results"

## Sample size used for creating density plots (evalMCMC.R)
densityPlotSampleSize <- 500

## Sample size used for creating density plots (predictMCMC.R)
quantilePlotSampleSize <- 50

hosp_nbinom_size = 35
FitTotalPeriod = 90

x11(width=20, height=12)
