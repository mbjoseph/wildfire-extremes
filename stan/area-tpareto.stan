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

  /**
  * Return the log probability of a unit-scale intrinsic
  * autoregressive (IAR) prior
  * with a sparse representation for the adjacency matrix
  *
  * @param phi Vector containing the parameters with a IAR prior
  * @param W_sparse Sparse representation of adjacency matrix (int array)
  * @param D_sparse Number of neighbors for each location (vector)
  * @param N Length of phi (int)
  * @param W_n Number of adjacent pairs (int)
  *
  * @return Log probability density of CAR prior up to additive constant
  */
  real sparse_iar_lpdf(vector phi,
    int[,] W_sparse, vector D_sparse, int N, int W_n) {
      row_vector[N] phit_D; // phi' * D
      row_vector[N] phit_W; // phi' * W

      phit_D = (phi .* D_sparse)';
      phit_W = rep_row_vector(0, N);
      for (i in 1:W_n) {
        phit_W[W_sparse[i, 1]] = phit_W[W_sparse[i, 1]] + phi[W_sparse[i, 2]];
        phit_W[W_sparse[i, 2]] = phit_W[W_sparse[i, 2]] + phi[W_sparse[i, 1]];
      }

      return 0.5 * -(phit_D * phi - (phit_W * phi));
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

  matrix<lower = 0, upper = 1>[N, N] W; // adjacency matrix
  int W_n;                // number of adjacent region pairs
  int<lower = 1, upper = N*T> tb_idx[n_u_tb - 1];
}

transformed data {
  int n_burn_mu = n_u_tb - 1;
  vector[n_fire] raw_sizes;
  vector[n_holdout_b] raw_holdout_b;
  int W_sparse[W_n, 2];   // adjacency pairs
  vector[N] D_sparse;     // diagonal of D (number of neigbors for each site)

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

  // iar params
  vector<lower = 0, upper = 1>[M] eta;   // autoregressive param
  vector<lower = 0>[M] sigma_phi;        // time difference spatial sd
  matrix[T, N] phiR;                     // unscaled values
}

transformed parameters {
  vector<lower = 0>[n_burn_mu] mu_burn[M];
  matrix[M, p] beta;
  matrix[M, p] lambda_tilde;
  vector[p] lambda_sq;
  vector[M] c;
  matrix[T, N] phi;
  vector[N * T] phi_vec;


  phi[1] = (phiR[1] - mean(phiR[1])) * sigma_phi[1];
  for (t in 2:T)
    phi[t] = eta[1] * phi[t -1] + (phiR[t] - mean(phiR[t])) * sigma_phi[1];

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

  for (i in 1:M)
    mu_burn[i] = exp(alpha[i]
               + csr_matrix_times_vector(n_burn_mu, p, w_tb, v_tb, u_tb, beta[i, ]')
               + phi_vec[tb_idx]);
}

model {
  lambda ~ cauchy(0, 1);
  L_beta ~ lkj_corr_cholesky(3);
  to_vector(betaR) ~ normal(0, 1);

  c_aux ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);
  alpha ~ normal(0, 5);
  tau ~ normal(0, 5);
  theta ~ cauchy(0, 1);

  eta ~ beta(8, 2);
  sigma_phi ~ normal(0, 1);
  for (t in 1:T) {
    phiR[t]' ~ sparse_iar(W_sparse, D_sparse, N, W_n);
  }

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
      mu_full[i] = exp(alpha[i]
                      + csr_matrix_times_vector(N * T, p, w, v, u, beta[i, ]')
                      + phi_vec);
  }

  for (i in 1:n_holdout_b) {
    holdout_loglik_b[i] = tpareto_lpdf(raw_holdout_b[i] | min_size,
                                                            mu_full[1][holdout_b_idx[i]],
                                                            theta);
    holdout_rep[i] = tpareto_rng(min_size, mu_full[1][holdout_b_idx[i]], theta) - min_size;
  }
}
