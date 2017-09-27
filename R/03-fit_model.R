stan_d <- list(
  N = N,
  T = T,
  L = 2,
  p = ncol(X),
  p_c = ncol(Xc),

  n_w = length(sparse_X$w),
  w = sparse_X$w,
  v = sparse_X$v,
  u = sparse_X$u,
  nrowX = nrow(X),
  n_u = length(sparse_X$u),
  burn_idx = burn_idx,

  n_wc = length(sparse_Xc$w),
  wc = sparse_Xc$w,
  vc = sparse_Xc$v,
  uc = sparse_Xc$u,

  counts = train_counts$n_fire,
  n_fire = nrow(train_burns),
  sizes = log(train_burns$R_ACRES - 1e3),
  log_area = log(train_counts$area * 1e-10)
)

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
                  iter = 1000)
# write_rds(m_fit, 'm_fit.rds')
