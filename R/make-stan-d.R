
# Bundle up the data for Stan models --------------------------------------


assert_that(identical(levels(area_df$NA_L3NAME),
                      levels(factor(st_covs$NA_L3NAME))))

scale_adj <- 1e4

stan_d <- list(
  N = N,
  T = T,
  p = length(colnamesX),

  n_count = nrow(train_counts),
  counts = train_counts$n_fire,

  log_area = log(area_df$area * 1e-11),
  er_idx_train = as.numeric(factor(train_counts$NA_L3NAME,
                                   levels = levels(area_df$NA_L3NAME))),
  er_idx_full = as.numeric(factor(st_covs$NA_L3NAME)),

  n_fire = nrow(train_burns),
  sizes = (train_burns$R_ACRES - 1e3) / scale_adj,
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

  M = 5,
  slab_df = 5,
  slab_scale = 2,

  eps_idx_train = eps_idx_train,
  eps_idx_future = eps_idx_future,
  size_threshold = 0,

  n_holdout_c = length(holdout_c_idx),
  holdout_c_idx = holdout_c_idx,
  holdout_c = holdout_counts$n_fire,

  n_holdout_b = length(holdout_b_idx),
  holdout_b_idx = holdout_b_idx,
  holdout_b = (holdout_burns$R_ACRES - 1e3) / scale_adj
  )
