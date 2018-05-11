library(readr)
library(rstan)
zi_d <- read_rds('data/processed/zi_d.rds')

zinb_init <- stan_model('stan/counts-zinb.stan')
zinb_full_fit <- sampling(zinb_init,
                          data = zi_d,
                          init_r = 0.01,
                          iter = 1000, 
                          pars =  c('beta', 'tau', 'alpha', 
                                    'c', 'mu_full', 'Rho_beta',
                                     'count_pred',
                                     'lambda_tilde',
                                     'holdout_loglik_c', 'train_loglik_c',
                                     'sigma_phi', 'phi', 'eta'),
                          cores = 4)
write_rds(zinb_full_fit, path = 'zinb_full_fit.rds')
