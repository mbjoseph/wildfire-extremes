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

  int<lower = 0> n_edges;
  int<lower = 1, upper = N> node1[n_edges];
  int<lower = 1, upper = N> node2[n_edges];
  int<lower = 1, upper = N*T> tb_idx[n_u_tb - 1]; // train burn spatial idx
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
  real<lower = 0> scale;

  // iar params
  vector<lower = 0, upper = 1>[M] eta;   // autoregressive param
  vector<lower = 0>[M] sigma_phi;        // time difference spatial sd
  vector[N] phiR[T];                     // unscaled values
}

transformed parameters {
  vector[n_burn_mu] mu_burn;
  matrix[M, p] beta;
  matrix[M, p] lambda_tilde;
  vector[p] lambda_sq;
  vector[M] c;
  matrix[T, N] phi;
  vector[N * T] phi_vec;

  phi[1] = phiR[1]' * sigma_phi[1];
  for (t in 2:T)
    phi[t] = eta[1] * phi[t - 1] + phiR[t]' * sigma_phi[1];

  phi_vec = to_vector(phi');


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

  mu_burn = alpha[1]
            + csr_matrix_times_vector(n_burn_mu, p, w_tb, v_tb, u_tb, beta[1,]')
            + phi_vec[tb_idx];
}

model {
  lambda ~ cauchy(0, 1);
  L_beta ~ lkj_corr_cholesky(3);
  to_vector(betaR) ~ normal(0, 1);

  c_aux ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);
  alpha ~ normal(0, 5);
  tau ~ normal(0, 5);
  scale ~ normal(0, 5);

  eta ~ beta(8, 2);
  sigma_phi ~ normal(0, 1);
  for (t in 1:T) {
    // IAR prior
    target += -0.5 * dot_self(phiR[t][node1] - phiR[t][node2]);
    sum(phiR[t]) ~ normal(0, .001 * N);
  }

  // fire sizes
  sizes ~ lognormal(mu_burn[burn_idx], scale);
}

generated quantities {
  vector[N * T] mu_full;
  vector[n_fire] loglik_f;
  vector[n_fire] size_rep;
  matrix[M, M] Rho_beta = multiply_lower_tri_self_transpose(L_beta);
  vector[n_holdout_b] holdout_loglik_b;
  vector[n_holdout_b] holdout_rep;

  for (i in 1:n_fire) {
    loglik_f[i] = lognormal_lpdf(sizes[i] | mu_burn[burn_idx[i]], scale);
    size_rep[i] = lognormal_rng(mu_burn[burn_idx[i]], scale);
  }

  // expected values
  mu_full = alpha[1] + csr_matrix_times_vector(N * T, p, w, v, u, beta[1, ]') + phi_vec;

  for (i in 1:n_holdout_b) {
    holdout_loglik_b[i] = lognormal_lpdf(holdout_b[i] | mu_full[holdout_b_idx[i]], scale);
    holdout_rep[i] = lognormal_rng(mu_full[holdout_b_idx[i]], scale);
  }
}
