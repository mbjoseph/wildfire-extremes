library(rstan)
library(shinystan)
options(mc.cores = 2)
source('R/eda.R')

stan_d <- list(n = nrow(extremes),
               y = extremes$max_y)

m_fit <- stan('stan/gev.stan', iter = 1000, data = stan_d, chains = 4,
              control = list(adapt_delta = .99))
traceplot(m_fit, inc_warmup = TRUE)
pairs(m_fit)

launch_shinystan(m_fit)
