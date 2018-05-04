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

  int<lower = 2, upper = 2> M; // num dimensions
  real<lower = 0> slab_df;
  real<lower = 0> slab_scale;

  // holdout counts
  int<lower = 1> n_holdout_c;
  int<lower = 1, upper = N * T> holdout_c_idx[n_holdout_c];
  int<lower = 0> holdout_c[n_holdout_c];

  int<lower = 1, upper = N * T> eps_idx_train[n_count];
  int<lower = 0> n_edges;
  int<lower = 1, upper = N> node1[n_edges];
  int<lower = 1, upper = N> node2[n_edges];
}

transformed data {
  vector[n_count] log_area_train = log_area[er_idx_train];
  vector[N * T] log_area_full = log_area[er_idx_full];
}

parameters {
  matrix[M, p] betaR;
  vector<lower = 0>[p] lambda;
  vector<lower = 0>[M] tau;
  vector[M] alpha; // intercept
  vector<lower = 0>[M] c_aux;
  cholesky_factor_corr[M] L_beta;
  real<lower = 0> nb_prec;

  // iar params
  vector<lower = 0, upper = 1>[M] eta;   // autoregressive param
  vector<lower = 0>[M] sigma_phi;        // time difference spatial sd
  vector[N] phiR[M, T];                  // unscaled values
}

transformed parameters {
  vector[n_count] mu_count;
  vector[n_count] logit_p;
  matrix[M, p] beta;
  matrix[M, p] lambda_tilde;
  vector[p] lambda_sq;
  vector[M] c;
  matrix[T, N] phi[M];
  vector[N * T] phi_vec[M];

  for (i in 1:M) {
    phi[i][1] = phiR[i, 1]' * sigma_phi[i];
    for (t in 2:T) {
      // subsequent timesteps
      phi[i][t] = eta[i] * phi[i][t - 1] + phiR[i, t]' * sigma_phi[i];
    }
    phi_vec[i] = to_vector(phi[i]');
  }

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

  mu_count = alpha[1]
             + csr_matrix_times_vector(n_count, p, w_tc, v_tc, u_tc, beta[1, ]')
             + phi_vec[1][eps_idx_train]
             + log_area_train;
  logit_p = alpha[2]
             + csr_matrix_times_vector(n_count, p, w_tc, v_tc, u_tc, beta[2, ]')
             + phi_vec[2][eps_idx_train];
}

model {
  lambda ~ cauchy(0, 1);
  L_beta ~ lkj_corr_cholesky(3);
  to_vector(betaR) ~ normal(0, 1);

  c_aux ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);
  alpha ~ normal(0, 5);
  tau ~ normal(0, 5);
  nb_prec ~ normal(0, 5);

  sigma_phi ~ normal(0, 1);
  eta ~ beta(8, 2);
  for (i in 1:M) {
    for (t in 1:T) {
      // IAR prior
      target += -0.5 * dot_self(phiR[i, t][node1] - phiR[i, t][node2]);
      sum(phiR[i, t]) ~ normal(0, .001 * N);
    }
  }

  // number of fires
  for (i in 1:n_count) {
    if (counts[i] == 0) {
      target += log_sum_exp(bernoulli_logit_lpmf(0 | logit_p[i]),
                            bernoulli_logit_lpmf(1 | logit_p[i])
                            + neg_binomial_2_log_lpmf(counts[i] | mu_count[i], nb_prec));
    } else {
      target += bernoulli_logit_lpmf(1 | logit_p[i])
                + neg_binomial_2_log_lpmf(counts[i] | mu_count[i], nb_prec);
    }
  }
}

generated quantities {
  vector[N * T] mu_full[M];
  matrix[M, M] Rho_beta = multiply_lower_tri_self_transpose(L_beta);
  vector[N * T] count_pred;
  vector[n_count] pearson_resid;
  vector[n_holdout_c] holdout_loglik_c;
  vector[n_count] train_loglik_c;

  // expected values
  for (i in 1:M){
    mu_full[i] = alpha[i]
               + csr_matrix_times_vector(N * T, p, w, v, u, beta[i, ]')
                + phi_vec[i];
  }

  // expected log counts need an offset for area
  mu_full[1] = mu_full[1] + log_area_full;

  // posterior predictions for the number of fires
  for (i in 1:(N * T)) {
    count_pred[i] = bernoulli_logit_rng(mu_full[2][i]) * neg_binomial_2_log_rng(mu_full[1][i], nb_prec);
  }

  {
    // compute pearson residuals
    vector[N * T] pr = inv_logit(mu_full[2]);
    vector[N * T] mu_nb = exp(mu_full[1]);
    vector[N * T] var_pr = pr .* (1 - pr);
    vector[N * T] var_nb = mu_nb + square(mu_nb) ./ nb_prec;
    vector[N * T] zinb_mean = pr .* mu_nb;
    vector[N * T] zinb_var = square(pr) .* var_nb + square(mu_nb) .* var_pr + var_nb .* var_pr;
    for (i in 1:n_count)
      pearson_resid[i] = (counts[i] - zinb_mean[eps_idx_train[i]]) / sqrt(zinb_var[eps_idx_train[i]]);
  }

  // training log likelihoods
  for (i in 1:n_count) {
    if (counts[i] == 0) {
      train_loglik_c[i] = log_sum_exp(bernoulli_logit_lpmf(0 | logit_p[i]),
                            bernoulli_logit_lpmf(1 | logit_p[i])
                            + neg_binomial_2_log_lpmf(counts[i] | mu_count[i], nb_prec));
    } else {
      train_loglik_c[i] = bernoulli_logit_lpmf(1 | logit_p[i])
                + neg_binomial_2_log_lpmf(counts[i] | mu_count[i], nb_prec);
    }
  }

  // holdout log likelihoods
  for (i in 1:n_holdout_c) {
    if (holdout_c[i] == 0) {
      holdout_loglik_c[i] = log_sum_exp(bernoulli_logit_lpmf(0 | mu_full[2][holdout_c_idx[i]]),
                            bernoulli_logit_lpmf(1 | mu_full[2][holdout_c_idx[i]])
                            + neg_binomial_2_log_lpmf(holdout_c[i] | mu_full[1][holdout_c_idx[i]], nb_prec));
    } else {
      holdout_loglik_c[i] = bernoulli_logit_lpmf(1 | mu_full[2][holdout_c_idx[i]])
                + neg_binomial_2_log_lpmf(holdout_c[i] | mu_full[1][holdout_c_idx[i]], nb_prec);
    }
  }
}
