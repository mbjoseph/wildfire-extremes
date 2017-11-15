source('R/02-explore.R')
source('R/make-stan-d.R')

pars <- c('beta', 'tau', 'alpha', 'c', 'mu_full', 'Rho_beta',
          'loglik_c', 'loglik_f')

control_list <- list(
  adapt_delta = 0.8,
  max_treedepth = 10
)

n_iter <- 1000

# Fancy models --------------------------------------------------------

w_init <- stan_model('stan/st-basis-nb-weibull.stan')
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


g_init <- stan_model('stan/st-basis-nb-gamma.stan')
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


ln_init <- stan_model('stan/st-basis-nb-lognorm.stan')
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

gpd_init <- stan_model('stan/st-basis-nb-gpd.stan')
gpd_fit <- sampling(
  gpd_init,
  data = gpd_d,
  cores = 4,
  init_r = 0.01,
  iter = n_iter,
  refresh = 1,
  control = control_list
)
write_rds(gpd_fit, paste0('gpdfit_',
                         Sys.time() %>% gsub(' ', '_', x = .),
                         '.rds'))
rm(gpd_fit)
gc()


zinb_gpd_init <- stan_model('stan/st-basis-zinb-gpd.stan')
zinb_gpd_fit <- sampling(
  zinb_gpd_init,
  data = zinb_gpd_d,
  cores = 4,
  init_r = 0.01,
  iter = n_iter,
  refresh = 1,
  control = control_list
)
write_rds(zinb_gpd_fit, paste0('zinbgpdfit_',
                          Sys.time() %>% gsub(' ', '_', x = .),
                          '.rds'))
