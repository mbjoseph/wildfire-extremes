functions {
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

  // holdout burn areas
  int<lower = 1> n_holdout_b;
  int<lower = 1, upper = N * T> holdout_b_idx[n_holdout_b];
  vector[n_holdout_b] holdout_b;

  matrix<lower = 0, upper = 1>[N, N] W; // adjacency matrix
  int W_n;                // number of adjacent region pairs
  int<lower = 1, upper = N*T> tb_idx[n_u_tb - 1];
}

transformed data {
  int n_burn_mu = n_u_tb - 1;
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
}


parameters {
  matrix[M, p] betaR;
  vector<lower = 0>[p] lambda;
  vector<lower = 0>[M] tau;
  vector[M] alpha; // intercept
  vector<lower = 0>[M] c_aux;
  cholesky_factor_corr[M] L_beta;
  real<lower = 0> shape;

  // iar params
  vector<lower = 0, upper = 1>[M] eta;   // autoregressive param
  vector<lower = 0>[M] sigma_phi;        // time difference spatial sd
  matrix[T, N] phiR;                     // unscaled values

}

transformed parameters {
  vector<lower = 0>[n_burn_mu] mu_burn;
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

  mu_burn = exp(alpha[1]
              + csr_matrix_times_vector(n_burn_mu, p, w_tb, v_tb, u_tb, beta[1,]')
              + phi_vec[tb_idx]);
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
  sizes ~ weibull(shape, mu_burn[burn_idx]);
}

generated quantities {
  vector[N * T] mu_full;
  vector[n_fire] loglik_f;
  vector[n_fire] size_rep;
  matrix[M, M] Rho_beta = multiply_lower_tri_self_transpose(L_beta);
  vector[n_holdout_b] holdout_loglik_b;
  vector[n_holdout_b] holdout_rep;

  for (i in 1:n_fire) {
    loglik_f[i] = weibull_lpdf(sizes[i] | shape, mu_burn[burn_idx[i]]);
    size_rep[i] = weibull_rng(shape, mu_burn[burn_idx[i]]);
  }

  mu_full = exp(alpha[1] + csr_matrix_times_vector(N * T, p, w, v, u, beta[1, ]') + phi_vec);
  for (i in 1:n_holdout_b) {
    holdout_loglik_b[i] = weibull_lpdf(holdout_b[i] | shape, mu_full[holdout_b_idx[i]]);
    holdout_rep[i] = weibull_rng(shape, mu_full[holdout_b_idx[i]]);
  }
}
