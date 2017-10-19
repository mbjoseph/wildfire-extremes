data {
  int<lower = 1> N; // # spatial units
  int<lower = 1> T; // # timesteps
  int<lower = 1> p; // # columns in design matrix

  int<lower = 1> n_count; // number of counts in training set
  int<lower = 0> counts[n_count]; // # of events in each spatial unit * timestep

  // areas of each ecoregion with indices for broadcasting
  vector[N] log_area;
  int<lower = 1, upper = N> er_idx_train[n_count];
  int<lower = 1, upper = N> er_idx_full[N * T];

  int n_fire;
  vector[n_fire] sizes;
  int<lower = 1, upper = N * T> burn_idx[n_fire];

  // full sparse matrix for all ecoregions X timesteps
  int<lower = 1> n_w;
  vector[n_w] w;
  int<lower = 1> v[n_w];
  int<lower = 1> u[N * T + 1];

  // training design matrix for counts
  int<lower = 1> n_w_tc;
  vector[n_w_tc] w_tc;
  int<lower = 1> v_tc[n_w_tc];
  int<lower = 1, upper = N * T + 1> n_u_tc;
  int<lower = 1> u_tc[n_u_tc];

  // training design matrix for burn areas
  int<lower = 1> n_w_tb;
  vector[n_w_tb] w_tb;
  int<lower = 1> v_tb[n_w_tb];
  int<lower = 1, upper = N * T + 1> n_u_tb;
  int<lower = 1> u_tb[n_u_tb];

  int<lower = 1, upper = 3> M; // num dimensions
  real<lower = 0> slab_df;
  real<lower = 0> slab_scale;
}

transformed data {
  vector[n_count] log_area_train = log_area[er_idx_train];
  vector[N * T] log_area_full = log_area[er_idx_full];
  int n_burn_mu = n_u_tb - 1;
}

parameters {
  matrix[M, p] betaR;
  vector<lower = 0>[p] lambda;
  vector<lower = 0>[M] tau;
  vector[M] alpha; // intercept
  vector<lower = 0>[M] c_aux;
  cholesky_factor_corr[M] L_beta;
  real<lower = 0> phi;
}

transformed parameters {
  vector[n_count] mu_count;
  vector[n_burn_mu] mu_burn[2];
  matrix[M, p] beta;
  matrix[M, p] lambda_tilde;
  vector[p] lambda_sq;
  vector<lower = 0>[M] c;

  c = slab_scale * sqrt(c_aux);

  // regularized horseshoe prior
  lambda_sq = square(lambda);
  for (i in 1:M) {
    lambda_tilde[i, ] = sqrt(
                            c[i]^2 * lambda_sq ./
                            (c[i]^2 + tau[i]^2 * lambda_sq)
                          )';
  }

  // multivariate horseshoe
  beta = diag_pre_multiply(tau, L_beta) *  betaR .* lambda_tilde;

  for (i in 1:2)
    mu_burn[i] = alpha[i] +
               + csr_matrix_times_vector(n_burn_mu, p, w_tb, v_tb, u_tb, beta[i, ]');

   mu_count = alpha[3]
             + csr_matrix_times_vector(n_count, p, w_tc, v_tc, u_tc, beta[3, ]')
             + log_area_train;
}

model {
  lambda ~ cauchy(0, 1);
  L_beta ~ lkj_corr_cholesky(2);
  to_vector(betaR) ~ normal(0, 1);

  c_aux ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);
  alpha ~ normal(0, 5);
  tau ~ normal(0, 1);
  phi ~ normal(0, 5);

  // number of fires
  counts ~ neg_binomial_2_log(mu_count, phi);

  // fire sizes
  sizes ~ weibull(exp(mu_burn[1])[burn_idx], exp(mu_burn[2])[burn_idx]);
}

generated quantities {
  vector[N * T] mu_full[M];
  vector[n_count] loglik_c;
  vector[n_fire] loglik_f;
  matrix[M, M] Rho_beta = multiply_lower_tri_self_transpose(L_beta);

  for (i in 1:n_count) {
    loglik_c[i] = neg_binomial_2_log_lpmf(counts[i] | mu_count[i], phi);
  }

  for (i in 1:n_fire) {
    loglik_f[i] = weibull_lpdf(sizes[i] | exp(mu_burn[1][burn_idx[i]]), exp(mu_burn[2][burn_idx[i]]));
  }

  // expected values
  for (i in 1:M){
    mu_full[i] = alpha[i]
               + csr_matrix_times_vector(N * T, p, w, v, u, beta[i, ]');

  }

  // expected log counts need an offset for area
  mu_full[M] = mu_full[M] + log_area_full;
}
