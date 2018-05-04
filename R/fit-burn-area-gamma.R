library(readr)
library(rstan)
stan_d <- read_rds('data/processed/stan_d.rds')

ba_gamma_init <- stan_model('stan/area-gamma.stan')
ba_gamma_fit <- sampling(
  ba_gamma_init,
  data = stan_d,
  cores = 4,
  init_r = .01,
  iter = 100)
write_rds(ba_gamma_fit, 'ba_gamma_fit.rds')
