library(readr)
library(rstan)
stan_d <- read_rds('data/processed/stan_d.rds')

ba_lognormal_init <- stan_model('stan/area-lognormal.stan')
ba_lognormal_fit <- sampling(
  ba_lognormal_init,
  data = stan_d,
  cores = 4,
  init_r = .01,
  iter = 1500)
write_rds(ba_lognormal_fit, 'ba_lognormal_fit.rds')
