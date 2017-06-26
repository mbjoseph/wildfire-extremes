data {
  int<lower = 1> n;
  int<lower = 1> T;
  int<lower = 1> L; // number of dimensions
  int<lower = 1> p; // number of columns in design matrix
  matrix[n * L, p] X[T];
  int<lower = 1> r;
  matrix[n * L, r] S[T];
  matrix[r, r] M[T];
  int<lower = 0> counts[n * T];
  int<lower = 1, upper = n * L * T> count_idx[n * T];
  int n_fire;
  vector[n_fire] sizes;
  int<lower = 1, upper = n * L * T> size_idx[n_fire];
}

parameters {
  vector[p] beta;
  matrix[r, T] etaR;
  cholesky_factor_corr[r] L_eta;
  real<lower = 0> sigma_0;
  vector<lower = 0>[r] sigma_eta;
  real<lower = 0> sigma_size;
  vector[n * L * T] epsilonR;
  real<lower = 0> sigma_mu;
}

transformed parameters {
  matrix[r, T] eta;
  vector[n * L * T] mu;

  eta[, 1] = diag_pre_multiply(rep_vector(sigma_0, r), L_eta) *  etaR[, 1];
  for (t in 2:T) {
    eta[, t] = M[t] * eta[, t - 1]
               + diag_pre_multiply(sigma_eta, L_eta) *  etaR[, t];
  }

  for (t in 1:T) {
    mu[(1 + (t - 1) * n * L):(t * n * L)] = X[t] * beta + S[t] * eta[, t];
  }
  // process error
  mu = mu + epsilonR * sigma_mu;
}

model {
  beta ~ normal(0, 5);
  to_vector(etaR) ~ normal(0, 1);
  L_eta ~ lkj_corr_cholesky(2);
  sigma_0 ~ normal(0, 3);
  sigma_eta ~ normal(0, 3);
  sigma_mu ~ normal(0, 3);
  epsilonR ~ normal(0, 1);

  // number of fires
  counts ~ poisson_log(mu[count_idx]);

  // fire sizes
  sigma_size ~ normal(0, 1);
  sizes ~ normal(mu[size_idx], sigma_size);
}

generated quantities {
  matrix[r, r] R_eta;

  R_eta = multiply_lower_tri_self_transpose(L_eta);
}
