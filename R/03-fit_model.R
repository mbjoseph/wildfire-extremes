source('R/02-explore.R')

assert_that(identical(levels(area_df$NA_L3NAME),
                      levels(factor(st_covs$NA_L3NAME))))

stan_d <- list(
  N = N,
  T = T,
  p = length(colnamesX),

  n_count = nrow(train_counts),
  counts = train_counts$n_fire,

  log_area = log(area_df$area * 1e-10),
  er_idx_train = as.numeric(factor(train_counts$NA_L3NAME,
                                   levels = levels(area_df$NA_L3NAME))),
  er_idx_full = as.numeric(factor(st_covs$NA_L3NAME)),

  n_fire = nrow(train_burns),
  sizes = log(train_burns$R_ACRES - 1e3),
  burn_idx = burn_idx,

  n_w = length(sparse_X$w),
  w = sparse_X$w,
  v = sparse_X$v,
  u = sparse_X$u,

  # sparse design matrix for training counts
  n_w_tc = length(sparse_X_tc$w),
  w_tc = sparse_X_tc$w,
  v_tc = sparse_X_tc$v,
  n_u_tc = length(sparse_X_tc$u),
  u_tc = sparse_X_tc$u,

  # sparse design matrix for training burns
  n_w_tb = length(sparse_X_tb$w),
  w_tb = sparse_X_tb$w,
  v_tb = sparse_X_tb$v,
  n_u_tb = length(sparse_X_tb$u),
  u_tb = sparse_X_tb$u,

  burn_eps_idx = burn_eps_idx,

  M = 3,
  slab_df = 2,
  slab_scale = 2,

  eps_idx_train = eps_idx_train,
  eps_idx_future = eps_idx_future)

pars <- c('beta', 'tau', 'alpha', 'sigma', 'c', 'mu_full', 'Rho_beta', 'Rho_eps',
          'loglik_c', 'loglik_f')


# Fancy model, nonstationary variance -------------------------------------
m_init <- stan_model('stan/st-basis-ns-lognorm.stan')
m_fit <- sampling(m_init,
                   data = stan_d,
                   pars = pars,
                   cores = 4,
                   init_r = 0.01,
                   iter = 1000,
                   refresh = 1)
write_rds(m_fit, 'm_fit.rds')

# Fancy Weibull model --------------------------------------------------------
weibull_scale_adj <- 1e4
stan_d$sizes <- train_burns$R_ACRES / weibull_scale_adj

w_init <- stan_model('stan/st-basis-ns-weibull.stan')
w_fit <- sampling(w_init,
                  data = stan_d,
                  pars = pars,
                  cores = 4,
                  init_r = 0.01,
                  iter = 800,
                  refresh = 1)
write_rds(w_fit, 'w_fit.rds')
