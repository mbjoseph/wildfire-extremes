source('R/02-explore.R')
library(rstan)

stan_d <- list(
  N = N,
  T = T,
  L = 2,
  p = ncol(X),
  p_c = ncol(Xc),
  log_area = log(burn_covs$area * 1e-10),

  n_count = nrow(train_counts),
  counts = train_counts$n_fire,
  count_idx = count_idx,

  n_fire = nrow(train_burns),
  sizes = log(train_burns$R_ACRES - 1e3),
  burn_idx = burn_idx,

  n_wc = length(sparse_Xc$w),
  wc = sparse_Xc$w,
  vc = sparse_Xc$v,
  uc = sparse_Xc$u,

  n_w = length(sparse_X$w),
  w = sparse_X$w,
  v = sparse_X$v,
  u = sparse_X$u)

m_init <- stan_model('stan/st-basis.stan')
m_fit <- sampling(m_init,
                  data = stan_d,
                  pars = c('beta', 'tau', 'sigma_size',
                           'mu', 'mu_counts',
                           'sigma_eps',
                           'beta_c', 'tau_c',
                           'alpha',
                           'c_sq'),
                  cores = 4,
                  init_r = .01,
                  iter = 1000,
                  refresh = 50)
# write_rds(m_fit, 'm_fit.rds')
