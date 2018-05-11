library(readr)
library(rstan)
stan_d <- read_rds('data/processed/stan_d.rds')

pois_init <- stan_model('stan/counts-pois.stan')
pois_fit <- vb(pois_init,
               data = stan_d,
               eta = .2,
               init = 0,
               tol_rel_obj = 0.008,
               adapt_engaged = FALSE)
write_rds(pois_fit, path = 'pois_fit.rds')
