source('R/02-explore.R')
source('R/make-stan-d.R')

count_pars <- c('beta', 'tau', 'alpha', 'c', 'mu_full', 'Rho_beta',
                'count_pred',
                'lambda_tilde',
                'holdout_loglik_c', 'train_loglik_c',
                'sigma_phi', 'phi', 'eta')

burn_area_pars <- c('beta', 'tau', 'alpha', 'c', 'mu_full', 'Rho_beta',
                    'loglik_f',
                    'lambda_tilde',
                    'holdout_loglik_b', 'size_rep', 'holdout_rep',
                    'sigma_phi', 'phi', 'eta')

control_list <- list(
  adapt_delta = 0.99,
  max_treedepth = 11
)

n_iter <- 1000
ref_rate <- 1

# Count models --------------------------------------------------------
zi_d <- stan_d
zi_d$M <- 2

pois_init <- stan_model('stan/counts-pois.stan')
pois_fit <- vb(pois_init,
               data = stan_d,
               eta = .2,
               init = 0,
               pars = count_pars,
               tol_rel_obj = 0.008,
               adapt_engaged = FALSE)
write_rds(pois_fit, path = 'pois_fit.rds')


zip_init <- stan_model('stan/counts-zip.stan')
zip_fit <- vb(zip_init,
              data = zi_d,
              eta = .5,
              pars = count_pars,
              tol_rel_obj = 0.008,
              init = 0,
              adapt_engaged = FALSE)
write_rds(zip_fit, path = 'zip_fit.rds')

nb_init <- stan_model('stan/counts-nb.stan')
nb_fit <- vb(
  nb_init,
  data = stan_d,
  eta = .1,
  pars = c(count_pars, 'nb_prec'),
  tol_rel_obj = 0.008,
  init = 0,
  adapt_engaged = FALSE)
write_rds(nb_fit, path = 'nb_fit.rds')

zinb_init <- stan_model('stan/counts-zinb.stan')
zinb_fit <- vb(zinb_init,
              data = zi_d,
              eta = .1,
              pars = c(count_pars, 'nb_prec'),
              tol_rel_obj = 0.008,
              init = 0,
              adapt_engaged = FALSE)
write_rds(zinb_fit, path = 'zinb_fit.rds')


# Burn area models --------------------------------------------------------
ba_gamma_init <- stan_model('stan/area-gamma.stan')
ba_gamma_fit <- sampling(
  ba_gamma_init,
  data = stan_d,
  pars = burn_area_pars,
  cores = 4,
  refresh = ref_rate,
  init_r = .01,
  control = control_list,
  iter = n_iter)
write_rds(ba_gamma_fit, 'ba_gamma_fit.rds')
rm(ba_gamma_fit)
gc()

ba_pareto_init <- stan_model('stan/area-pareto.stan')
ba_pareto_fit <- sampling(
  ba_pareto_init,
  pars = c('shape', burn_area_pars),
  data = stan_d,
  cores = 4,
  refresh = ref_rate,
  init_r = .1,
  control = control_list,
  iter = n_iter)
write_rds(ba_pareto_fit, 'ba_pareto_fit.rds')

ba_tpareto_init <- stan_model('stan/area-tpareto.stan')
ba_tpareto_fit <- sampling(
  ba_tpareto_init,
  pars = burn_area_pars,
  data = stan_d,
  cores = 4,
  refresh = ref_rate,
  init_r = .01,
  control = control_list,
  iter = n_iter)
write_rds(ba_tpareto_fit, 'ba_tpareto_fit.rds')


ba_lognormal_init <- stan_model('stan/area-lognormal.stan')
ba_lognormal_fit <- sampling(
  ba_lognormal_init,
  pars = burn_area_pars,
  data = stan_d,
  cores = 4,
  refresh = ref_rate,
  init_r = .01,
  control = control_list,
  iter = n_iter)
write_rds(ba_lognormal_fit, 'ba_lognormal_fit.rds')


ba_weibull_init <- stan_model('stan/area-weibull.stan')
ba_weibull_fit <- sampling(
  ba_weibull_init,
  data = stan_d,
  pars = burn_area_pars,
  cores = 4,
  refresh = ref_rate,
  init_r = .01,
  control = control_list,
  iter = n_iter)
write_rds(ba_weibull_fit, 'ba_weibull_fit.rds')

