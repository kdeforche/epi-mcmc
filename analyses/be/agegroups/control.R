settings <- "../settings.R"

if (exists("total_deaths_at_lockdown")) {
  init <- c(3.6, 0.7, 3.04, 0.12, 0.27, 0.09,
            log(0.01), log(0.3), 14, 13, 14, 13, 1.3, 7, total_deaths_at_lockdown, 0)
  scales <- c(0.15, 0.05, 0.15, 0.05, 0.15, 0.05,
              0.05, 0.1, 1, 1, 1, 1, 0.3, 0.3, total_deaths_at_lockdown / 20, 0.1)
  scale <- 0.045 * scales
}

hosp_nbinom_size = 35
#mort_nbinom_size = 20
FitTotalPeriod <- 90
