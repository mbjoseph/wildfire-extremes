library(readr)
library(rstan)
zi_d <- read_rds('data/processed/zi_d.rds')

zip_init <- stan_model('stan/counts-zip.stan')
zip_fit <- vb(zip_init,
              data = zi_d,
              eta = .4,
              tol_rel_obj = 0.008,
              init = 0,
              adapt_engaged = FALSE)
write_rds(zip_fit, path = 'zip_fit.rds')
