library(readr)
library(rstan)
stan_d <- read_rds('data/processed/stan_d.rds')

ba_pareto_init <- stan_model('stan/area-pareto.stan')
ba_pareto_fit <- sampling(
  ba_pareto_init,
  data = stan_d,
  cores = 4,
  init_r = .01,
  iter = 1500)
write_rds(ba_pareto_fit, 'ba_pareto_fit.rds')
