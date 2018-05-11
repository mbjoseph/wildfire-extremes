library(readr)
library(rstan)
zi_d <- read_rds('data/processed/zi_d.rds')

zinb_init <- stan_model('stan/counts-zinb.stan')
zinb_fit <- vb(zinb_init,
               data = zi_d,
               eta = .3,
               tol_rel_obj = 0.008,
               init = 0,
               adapt_engaged = FALSE)
write_rds(zinb_fit, path = 'zinb_fit.rds')
