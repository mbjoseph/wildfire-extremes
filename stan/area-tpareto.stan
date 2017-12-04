functions {
  real tpareto_lpdf(real y, real a, real lambda, real theta) {
    // tapered Pareto log pdf
      return log(lambda / y + 1 / theta) + lambda * (log(a) - log(y)) + (a - y) / theta;
  }
  real tpareto_rng(real a, real lambda, real theta) {
    {
      real x;
      real w;
      // tapered Pareto rng
      x = pareto_rng(a, lambda);
      w = exponential_rng(1 / theta) + a;
      if (x < w) {
        return x;
      } else {
        return w;
      }
    }
  }
}
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

  real<lower = 0> size_threshold;

  // holdout burn areas
  int<lower = 1> n_holdout_b;
  int<lower = 1, upper = N * T> holdout_b_idx[n_holdout_b];
  vector[n_holdout_b] holdout_b;

  real<lower = 0> min_size;
}

transformed data {
  int n_burn_mu = n_u_tb - 1;
  vector[n_fire] raw_sizes;
  vector[n_holdout_b] raw_holdout_b;

  raw_sizes = sizes + min_size;
  raw_holdout_b = holdout_b + min_size;
}

parameters {
  matrix[M, p] betaR;
  vector<lower = 0>[p] lambda;
  vector<lower = 0>[M] tau;
  vector[M] alpha; // intercept
  vector<lower = 0>[M] c_aux;
  cholesky_factor_corr[M] L_beta;
  real<lower = 0> theta;
}

transformed parameters {
  vector<lower = 0>[n_burn_mu] mu_burn[M];
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

  for (i in 1:M)
    mu_burn[i] = exp(alpha[i]
               + csr_matrix_times_vector(n_burn_mu, p, w_tb, v_tb, u_tb, beta[i, ]'));
}

model {
  lambda ~ cauchy(0, 1);
  L_beta ~ lkj_corr_cholesky(3);
  to_vector(betaR) ~ normal(0, 1);

  c_aux ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);
  alpha ~ normal(0, 5);
  tau ~ normal(0, 5);
  theta ~ cauchy(0, 1);

  // fire sizes
  for (i in 1:n_fire)
    raw_sizes[i] ~ tpareto(min_size, mu_burn[1][burn_idx[i]], theta);
}

generated quantities {
  vector<lower = 0>[N * T] mu_full[M];
  vector[n_fire] loglik_f;
  vector[n_fire] size_rep;
  matrix[M, M] Rho_beta = multiply_lower_tri_self_transpose(L_beta);
  vector[n_holdout_b] holdout_loglik_b;
  vector[n_holdout_b] holdout_rep;

  for (i in 1:n_fire) {
    loglik_f[i] = tpareto_lpdf(raw_sizes[i] | min_size,
                                                mu_burn[1][burn_idx[i]],
                                                theta);
    size_rep[i] = tpareto_rng(min_size, mu_burn[1][burn_idx[i]], theta) - min_size;
  }

  // expected values
  for (i in 1:M) {
      mu_full[i] = exp(alpha[i] + csr_matrix_times_vector(N * T, p, w, v, u, beta[i, ]'));
  }

  for (i in 1:n_holdout_b) {
    holdout_loglik_b[i] = tpareto_lpdf(raw_holdout_b[i] | min_size,
                                                            mu_full[1][holdout_b_idx[i]],
                                                            theta);
    holdout_rep[i] = tpareto_rng(min_size, mu_full[1][holdout_b_idx[i]], theta) - min_size;
  }
}
