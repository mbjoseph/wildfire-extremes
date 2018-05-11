library(readr)
library(rstan)
stan_d <- read_rds('data/processed/stan_d.rds')

ba_tpareto_init <- stan_model('stan/area-tpareto.stan')
ba_tpareto_fit <- sampling(
  ba_tpareto_init,
  data = stan_d,
  cores = 4,
  init_r = .01,
  iter = 1000)
write_rds(ba_tpareto_fit, 'ba_tpareto_fit.rds')
