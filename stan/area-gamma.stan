data {
  int<lower = 1> N; // # spatial units
  int<lower = 1> T; // # timesteps
  int<lower = 1> p; // # columns in design matrix

  int n_fire;
  vector[n_fire] sizes;
  int<lower = 1, upper = N * T> burn_idx[n_fire];

  // full sparse matrix for all ecoregions X timesteps
  int<lower = 1> n_w;
  vector[n_w] w;
  int<lower = 1> v[n_w];
  int<lower = 1> u[N * T + 1];
  // training design matrix for burn areas
  int<lower = 1> n_w_tb;
  vector[n_w_tb] w_tb;
  int<lower = 1> v_tb[n_w_tb];
  int<lower = 1, upper = N * T + 1> n_u_tb;
  int<lower = 1> u_tb[n_u_tb];

  int<lower = 1, upper = 1> M; // num dimensions
  real<lower = 0> slab_df;
  real<lower = 0> slab_scale;

  // holdout burn areas
  int<lower = 1> n_holdout_b;
  int<lower = 1, upper = N * T> holdout_b_idx[n_holdout_b];
  vector[n_holdout_b] holdout_b;
}

transformed data {
  int n_burn_mu = n_u_tb - 1;
}

parameters {
  matrix[M, p] betaR;
  vector<lower = 0>[p] lambda;
  vector<lower = 0>[M] tau;
  vector[M] alpha; // intercept
  vector<lower = 0>[M] c_aux;
  cholesky_factor_corr[M] L_beta;
  real<lower = 0> shape;
}

transformed parameters {
  vector<lower = 0>[n_burn_mu] mu_burn;
  matrix[M, p] beta;
  matrix[M, p] lambda_tilde;
  vector[p] lambda_sq;
  vector[M] c;

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

  mu_burn = exp(alpha[1] + csr_matrix_times_vector(n_burn_mu, p, w_tb, v_tb, u_tb, beta[1,]'));
}

model {
  lambda ~ cauchy(0, 1);
  L_beta ~ lkj_corr_cholesky(3);
  to_vector(betaR) ~ normal(0, 1);

  c_aux ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);
  alpha ~ normal(0, 5);
  tau ~ normal(0, 5);
  shape ~ normal(0, 5);

  // fire sizes
  sizes ~ gamma(shape, shape ./ mu_burn[burn_idx]);
}

generated quantities {
  vector<lower = 0>[N * T] mu_full;
  vector[n_fire] loglik_f;
  vector[n_fire] size_rep;
  matrix[M, M] Rho_beta = multiply_lower_tri_self_transpose(L_beta);
  vector[n_holdout_b] holdout_loglik_b;
  vector[n_holdout_b] holdout_rep;


  for (i in 1:n_fire) {
    loglik_f[i] = gamma_lpdf(sizes[i] | shape, shape ./ mu_burn[burn_idx[i]]);
    size_rep[i] = gamma_rng(shape, shape ./ mu_burn[burn_idx[i]]);
  }

  // expected values
  mu_full = exp(alpha[1] + csr_matrix_times_vector(N * T, p, w, v, u, beta[1, ]'));

  for (i in 1:n_holdout_b) {
    holdout_loglik_b[i] = gamma_lpdf(holdout_b[i] | shape,
                                                    shape ./ mu_full[holdout_b_idx[i]]);
    holdout_rep[i] = gamma_rng(shape, shape ./ mu_full[holdout_b_idx[i]]);
  }

}
