source('R/02-explore.R')

stan_d <- list(
  N = N,
  T = T,
  p = ncol(X),
  log_area = log(burn_covs$area * 1e-10),

  n_count = nrow(train_counts),
  counts = train_counts$n_fire,
  count_idx = count_idx,

  n_fire = nrow(train_burns),
  sizes = log(train_burns$R_ACRES - 1e3),
  burn_idx = burn_idx,

  n_w = length(sparse_X$w),
  w = sparse_X$w,
  v = sparse_X$v,
  u = sparse_X$u,
  M = 3,
  slab_df = 3,
  slab_scale = 5)



# Fancy model, nonstationary variance -------------------------------------
pars <- c('beta', 'tau', 'alpha', 'sigma', 'c', 'mu', 'Rho_beta', 'Rho_eps')
m_init <- stan_model('stan/st-basis-ns-lognorm.stan')
m_fit <- sampling(m_init,
                   data = stan_d,
                   pars = pars,
                   cores = 4,
                   init_r = 0.01,
                   iter = 800,
                   refresh = 1)
write_rds(m_fit, 'm_fit.rds')

# Fancy Weibull model --------------------------------------------------------
w_d <- list(
  N = N,
  T = T,
  L = 2,
  p = ncol(X),
  log_area = log(burn_covs$area * 1e-10),

  n_count = nrow(train_counts),
  counts = train_counts$n_fire,
  count_idx = count_idx,

  n_fire = nrow(train_burns),
  sizes = (train_burns$R_ACRES - 1e3) / 1e5,
  burn_idx = burn_idx,

  n_w = length(sparse_X$w),
  w = sparse_X$w,
  v = sparse_X$v,
  u = sparse_X$u)


w_init <- stan_model('stan/st-basis-ns-weibull.stan')
w_fit <- sampling(w_init,
                  data = w_d,
                  pars = pars,
                  cores = 4,
                  init_r = 0.01,
                  iter = 800,
                  refresh = 1)
write_rds(w_fit, 'w_fit.rds')
