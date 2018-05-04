library(readr)
library(rstan)
stan_d <- read_rds('data/processed/stan_d.rds')

ba_weibull_init <- stan_model('stan/area-weibull.stan')
ba_weibull_fit <- sampling(
  ba_weibull_init,
  data = stan_d,
  cores = 4,
  init_r = .01,
  iter = 100)
write_rds(ba_weibull_fit, 'ba_weibull_fit.rds')
