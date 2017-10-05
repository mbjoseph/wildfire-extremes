data {
  int<lower = 1> N; // # spatial units
  int<lower = 1> T; // # timesteps
  int<lower = 1> p; // # columns in design matrix
  vector[N * T] log_area;

  int<lower = 1> n_count; // number of counts in training set
  int<lower = 0> counts[n_count]; // # of events in each spatial unit * timestep
  int<lower = 1, upper = N * T> count_idx[n_count];

  int n_fire;
  vector[n_fire] sizes;
  int<lower = 1, upper = N * T> burn_idx[n_fire];

  // sparse matrix for fire sizes
  int<lower = 1> n_w;
  vector[n_w] w;
  int<lower = 1> v[n_w];
  int<lower = 1> u[N * T + 1];

  int<lower = 1, upper = 3> M; // num dimensions
  real<lower = 0> slab_df;
  real<lower = 0> slab_scale;
}

parameters {
  matrix[M, p] betaR;
  vector<lower = 0>[p] lambda[M];
  vector<lower = 0>[M] tau;
  vector[M] alpha; // intercept
  matrix[M, N * T] epsR;
  cholesky_factor_corr[M] L_eps;
  vector<lower = 0>[M] sigma;
  vector<lower = 0>[M] c_aux;
  cholesky_factor_corr[M] L_beta;
}

transformed parameters {
  vector[N * T] mu[M];
  matrix[M, p] beta;
  matrix[M, p] lambda_tilde;
  matrix[M, p] lambda_sq;
  vector<lower = 0>[M] c;
  matrix[M, N * T] eps;

  eps = diag_pre_multiply(sigma, L_eps) * epsR;

  c = slab_scale * sqrt(c_aux);

  // regularized horseshoe prior
  for (i in 1:M) {
    lambda_sq[i, ] = square(lambda[i])';
    lambda_tilde[i, ] = sqrt(
                            c[i]^2 * lambda_sq[i, ] ./
                            (c[i]^2 + tau[i]^2 * lambda_sq[i, ])
                          );
  }

  // multivariate horseshoe
  beta = diag_pre_multiply(tau, L_beta) *  betaR .* lambda_tilde;

  for (i in 1:M)
    mu[i] = alpha[i] +
              csr_matrix_times_vector(N * T, p, w, v, u, beta[i, ]') +
              eps[i, ]';

  // expected log counts need an offset for area
  mu[M] = mu[M] + log_area;
}

model {
  for (i in 1:M) {
    lambda[i] ~ cauchy(0, 1);
    epsR[i] ~ normal(0, 1);
  }
  L_beta ~ lkj_corr_cholesky(2);
  to_vector(betaR) ~ normal(0, 1);

  c_aux ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);
  sigma ~ normal(0, 1);
  L_eps ~ lkj_corr_cholesky(2);
  alpha ~ normal(0, 1);
  tau ~ normal(0, 1);

  // number of fires
  counts ~ poisson_log(mu[M][count_idx]);

  // fire sizes
  sizes ~ normal(mu[1][burn_idx], exp(mu[2])[burn_idx]);
}

generated quantities {
  matrix[M, M] Rho_beta = multiply_lower_tri_self_transpose(L_beta);
  matrix[M, M] Rho_eps = multiply_lower_tri_self_transpose(L_eps);
}
