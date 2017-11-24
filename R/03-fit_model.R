source('R/02-explore.R')
source('R/make-stan-d.R')

pars <- c('beta', 'tau', 'alpha', 'c', 'mu_full', 'Rho_beta',
          'loglik_c', 'loglik_f', 'count_pred',
          'lambda_tilde',
          'holdout_loglik_c',  'holdout_loglik_b')

control_list <- list(
  adapt_delta = 0.9,
  max_treedepth = 11
)

n_iter <- 1000

# Fancy models --------------------------------------------------------
zinb_lomax_init <- stan_model('stan/zinb-lomax.stan')

zinb_lomax_fit <- sampling(
  zinb_lomax_init,
  data = stan_d,
  cores = 4,
  init_r = 0.01,
  iter = n_iter,
  refresh = 1,
  control = control_list,
  pars = pars
)
write_rds(zinb_lomax_fit, paste0('zinblomaxfit_',
                               Sys.time() %>% gsub(' ', '_', x = .),
                               '.rds'))



w_init <- stan_model('stan/zinb-weibull.stan')
w_fit <- sampling(w_init,
                  data = stan_d,
                  pars = pars,
                  cores = 4,
                  init_r = 0.01,
                  iter = n_iter,
                  refresh = 1,
                  control = control_list)
write_rds(w_fit, paste0('wfit_',
                        Sys.time() %>% gsub(' ', '_', x = .),
                        '.rds'))
rm(w_fit)
gc()


g_init <- stan_model('stan/zinb-gamma.stan')
g_fit <- sampling(g_init,
                  data = stan_d,
                  pars = pars,
                  cores = 4,
                  init_r = 0.01,
                  iter = n_iter,
                  refresh = 1,
                  control = control_list)
write_rds(g_fit, paste0('gfit_',
                        Sys.time() %>% gsub(' ', '_', x = .),
                        '.rds'))
rm(g_fit)
gc()


ln_init <- stan_model('stan/zinb-lognorm.stan')
ln_fit <- sampling(ln_init,
                   data = stan_d,
                   pars = pars,
                   cores = 4,
                   init_r = 0.01,
                   iter = n_iter,
                   refresh = 1,
                   control = control_list)
write_rds(ln_fit, paste0('lnfit_',
                         Sys.time() %>% gsub(' ', '_', x = .),
                         '.rds'))
rm(ln_fit)
gc()

