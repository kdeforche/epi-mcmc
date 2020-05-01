source("model.R")

transformParams <- function(params)
{
    result = params
    result[4] = exp(params[3] + params[4])
    result[3] = exp(params[3])

    result
}

fit.paramnames <- c(fit.paramnames, "Nef")
keyparamnames <- c(keyparamnames, "Nef")
fitkeyparamnames <- c(fitkeyparamnames, "Nef")

init <- c(init, 0.8)
scales <- c(scales, 0.01)
