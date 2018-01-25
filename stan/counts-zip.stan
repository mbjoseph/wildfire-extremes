functions {
  /**
  * Return the log probability of a unit-scale proper
  * conditional autoregressive (CAR) prior
  * with a sparse representation for the adjacency matrix
  *
  * @param phi Vector containing the parameters with a CAR prior
  * @param gamma Dependence (usually spatial) parameter for the CAR prior (real)
  * @param W_sparse Sparse representation of adjacency matrix (int array)
  * @param N Length of phi (int)
  * @param W_n Number of adjacent pairs (int)
  * @param D_sparse Number of neighbors for each location (vector)
  * @param lambda_car Eigenvalues of D^{-1/2}*W*D^{-1/2} (vector)
  *
  * @return Log probability density of CAR prior up to additive constant
  */
  real sparse_car_lpdf(vector phi, real gamma,
    int[,] W_sparse, vector D_sparse, vector lambda_car, int N, int W_n) {
      row_vector[N] phit_D; // phi' * D
      row_vector[N] phit_W; // phi' * W
      vector[N] ldet_terms;

      phit_D = (phi .* D_sparse)';
      phit_W = rep_row_vector(0, N);
      for (i in 1:W_n) {
        phit_W[W_sparse[i, 1]] = phit_W[W_sparse[i, 1]] + phi[W_sparse[i, 2]];
        phit_W[W_sparse[i, 2]] = phit_W[W_sparse[i, 2]] + phi[W_sparse[i, 1]];
      }

      for (i in 1:N) ldet_terms[i] = log1m(gamma * lambda_car[i]);
      return 0.5 * (sum(ldet_terms) - (phit_D * phi - gamma * (phit_W * phi)));
  }
}
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
  matrix<lower = 0, upper = 1>[N, N] W; // adjacency matrix
  int W_n;                // number of adjacent region pairs
}

transformed data {
  vector[n_count] log_area_train = log_area[er_idx_train];
  vector[N * T] log_area_full = log_area[er_idx_full];
  int W_sparse[W_n, 2];   // adjacency pairs
  vector[N] D_sparse;     // diagonal of D (number of neigbors for each site)
  vector[N] lambda_car;       // eigenvalues of invsqrtD * W * invsqrtD

  { // generate sparse representation for W
  int counter;
  counter = 1;
  // loop over upper triangular part of W to identify neighbor pairs
    for (i in 1:(N - 1)) {
      for (j in (i + 1):N) {
        if (W[i, j] == 1) {
          W_sparse[counter, 1] = i;
          W_sparse[counter, 2] = j;
          counter = counter + 1;
        }
      }
    }
  }
  for (i in 1:N) D_sparse[i] = sum(W[i]);
  {
    vector[N] invsqrtD;
    for (i in 1:N) {
      invsqrtD[i] = 1 / sqrt(D_sparse[i]);
    }
    lambda_car = eigenvalues_sym(quad_form(W, diag_matrix(invsqrtD)));
  }
}

parameters {
  matrix[M, p] betaR;
  vector<lower = 0>[p] lambda;
  vector<lower = 0>[M] tau;
  vector[M] alpha; // intercept
  vector<lower = 0>[M] c_aux;
  cholesky_factor_corr[M] L_beta;

  // car params
  vector<lower = 0, upper = 1>[M] gamma; // spatial dependence
  vector<lower = 0, upper = 1>[M] eta;   // autoregressive param
  vector<lower = 0>[M] sigma_phi;        // time difference spatial sd
  matrix[N, T] phiR[M];                  // unscaled values
}


transformed parameters {
  vector[n_count] mu_count;
  vector[n_count] logit_p;
  matrix[M, p] beta;
  matrix[M, p] lambda_tilde;
  vector[p] lambda_sq;
  vector[M] c;
  matrix[N, T] phi[M];
  vector[N * T] phi_vec[M];

  for (i in 1:M) {
    for (j in 1:N) {
      phi[i][j, 1] = phiR[i][j, 1] * sigma_phi[i]; // initial time step
      for (t in 2:T) {
        // autoregressive structure for subsequent time steps
        phi[i][j, t] = eta[i] * phi[i][j, t - 1] + phiR[i][j, t] * sigma_phi[i];
      }
    }
    phi_vec[i] = to_vector(phi[i]);
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

  sigma_phi ~ normal(0, 1);
  for (i in 1:M) {
    for (t in 1:T) {
      phiR[i][, t] ~ sparse_car(gamma[i], W_sparse, D_sparse, lambda_car, N, W_n);
    }
  }

  // number of fires
  for (i in 1:n_count) {
    if (counts[i] == 0) {
      target += log_sum_exp(bernoulli_logit_lpmf(0 | logit_p[i]),
                            bernoulli_logit_lpmf(1 | logit_p[i])
                            + poisson_log_lpmf(counts[i] | mu_count[i]));
    } else {
      target += bernoulli_logit_lpmf(1 | logit_p[i])
                + poisson_log_lpmf(counts[i] | mu_count[i]);
    }
  }
}

generated quantities {
  vector[N * T] mu_full[M];
  matrix[M, M] Rho_beta = multiply_lower_tri_self_transpose(L_beta);
  vector[N * T] count_pred;
  vector[n_holdout_c] holdout_loglik_c;
  vector[n_count] train_loglik_c;
  vector[n_count] pearson_resid;

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
    count_pred[i] = bernoulli_logit_rng(mu_full[2][i]) * poisson_log_rng(mu_full[1][i]);
  }

    {
    // compute pearson residuals
    vector[N * T] pr = inv_logit(mu_full[2]);
    vector[N * T] mu_pois = exp(mu_full[1]);
    vector[N * T] var_pr = pr .* (1 - pr);
    vector[N * T] zip_mean = pr .* mu_pois;
    vector[N * T] zip_var = square(pr) .* mu_pois + square(mu_pois) .* var_pr + mu_pois .* var_pr;
    for (i in 1:n_count)
      pearson_resid[i] = (counts[i] - zip_mean[eps_idx_train[i]]) / sqrt(zip_var[eps_idx_train[i]]);
  }


  // training log likelihoods
  for (i in 1:n_count) {
    if (counts[i] == 0) {
      train_loglik_c[i] = log_sum_exp(bernoulli_logit_lpmf(0 | logit_p[i]),
                            bernoulli_logit_lpmf(1 | logit_p[i])
                            + poisson_log_lpmf(counts[i] | mu_count[i]));
    } else {
      train_loglik_c[i] = bernoulli_logit_lpmf(1 | logit_p[i])
                + poisson_log_lpmf(counts[i] | mu_count[i]);
    }
  }

  // holdout log likelihoods
  for (i in 1:n_holdout_c) {
    if (holdout_c[i] == 0) {
      holdout_loglik_c[i] = log_sum_exp(bernoulli_logit_lpmf(0 | mu_full[2][holdout_c_idx[i]]),
                            bernoulli_logit_lpmf(1 | mu_full[2][holdout_c_idx[i]])
                            + poisson_log_lpmf(counts[i] | mu_full[1][holdout_c_idx[i]]));
    } else {
      holdout_loglik_c[i] = bernoulli_logit_lpmf(1 | mu_full[2][holdout_c_idx[i]])
                + poisson_log_lpmf(counts[i] | mu_full[1][holdout_c_idx[i]]);
    }
  }
}
