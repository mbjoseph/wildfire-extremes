functions {
  /**
  * Return the log probability of a proper intrinsic autoregressive (IAR) prior
  * with a sparse representation for the adjacency matrix
  *
  * @param phi Vector containing the parameters with a IAR prior
  * @param tau Precision parameter for the IAR prior (real)
  * @param W_sparse Sparse representation of adjacency matrix (int array)
  * @param n Length of phi (int)
  * @param W_n Number of adjacent pairs (int)
  * @param D_sparse Number of neighbors for each location (vector)
  *
  * @return Log probability density of IAR prior up to additive constant
  */
  real sparse_iar_lpdf(vector phi,
    int[,] W_sparse, vector D_sparse, int n, int W_n) {
      row_vector[n] phit_D; // phi' * D
      row_vector[n] phit_W; // phi' * W
      vector[n] ldet_terms;

      phit_D = (phi .* D_sparse)';
      phit_W = rep_row_vector(0, n);
      for (i in 1:W_n) {
        phit_W[W_sparse[i, 1]] = phit_W[W_sparse[i, 1]] + phi[W_sparse[i, 2]];
        phit_W[W_sparse[i, 2]] = phit_W[W_sparse[i, 2]] + phi[W_sparse[i, 1]];
      }

      return 0.5 * (- (phit_D * phi - (phit_W * phi)));
  }
}

data {
  int<lower = 1> n;                       // sample size
  int<lower = 0> y[n];                    // number of fires
  int<lower = 1> J;                       // number of ecoregions
  int<lower = 1, upper = J> reg[n];       // ecoregion index
  vector[n] t;                            // time
  vector[n] log_offset;                   // log area
  matrix<lower = 0, upper = 1>[J, J] W;   // adjacency matrix
  int W_n;                                // number of adjacent region pairs
  int n_year;
  int year[n];
}

transformed data {
  int W_sparse[W_n, 2];   // adjacency pairs
  vector[J] D_sparse;     // diagonal of D (number of neigbors for each site)
  vector[J] lambda;       // eigenvalues of invsqrtD * W * invsqrtD

  { // generate sparse representation for W
  int counter;
  counter = 1;
  // loop over upper triangular part of W to identify neighbor pairs
    for (i in 1:(J - 1)) {
      for (j in (i + 1):J) {
        if (W[i, j] == 1) {
          W_sparse[counter, 1] = i;
          W_sparse[counter, 2] = j;
          counter = counter + 1;
        }
      }
    }
  }
  for (i in 1:J) D_sparse[i] = sum(W[i]);
  {
    vector[J] invsqrtD;
    for (i in 1:J) {
      invsqrtD[i] = 1 / sqrt(D_sparse[i]);
    }
    lambda = eigenvalues_sym(quad_form(W, diag_matrix(invsqrtD)));
  }
}

parameters {
  real alpha0;
  vector[J] phiR[n_year];
  vector<lower = 0>[n_year] sigma_phi;
  real mu_sd;
  real<lower = 0> sd_sd;
}

transformed parameters {
  matrix[J, n_year] phi;
  vector[n] log_lambda;

  phi[, 1] = phiR[1] * sigma_phi[1] + alpha0;
  for (i in 2:n_year) {
    phi[, i] = phi[, i - 1] + sigma_phi[i] * phiR[i];
  }

  log_lambda = log_offset;
  for (i in 1:n) {
    log_lambda[i] = log_lambda[i] + phi[reg[i], year[i]];
  }
}

model {
  for (i in 1:n_year) {
    phiR[i] ~ sparse_iar(W_sparse, D_sparse, J, W_n);
  }
  alpha0 ~ normal(0, 1);

  sigma_phi ~  lognormal(mu_sd, sd_sd);
  mu_sd ~ normal(0, 1);
  sd_sd ~ lognormal(0, 1);

  y ~ poisson_log(log_lambda);
}
