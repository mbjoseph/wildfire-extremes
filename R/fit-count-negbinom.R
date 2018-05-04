library(readr)
library(rstan)
stan_d <- read_rds('data/processed/stan_d.rds')

nb_init <- stan_model('stan/counts-nb.stan')
nb_fit <- vb(
  nb_init,
  data = stan_d,
  eta = .4,
  tol_rel_obj = 0.008,
  init = 0,
  adapt_engaged = FALSE)
write_rds(nb_fit, path = 'nb_fit.rds')
