source("model.R")

transformParams <- function(params)
{
    result = params
    result[4] = exp(params[3] + params[4])
    result[3] = exp(params[3])

    result
}

fit.paramnames <- append(fit.paramnames, "Nef", after=9)
keyparamnames <- append(keyparamnames, "Nef", after=9)
fitkeyparamnames <- append(fitkeyparamnames, "Nef", after=9)

init <- append(init, 0.8, after=9)
scales <- append(scales, 0.01, after=9)
