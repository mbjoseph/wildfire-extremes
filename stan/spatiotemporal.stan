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
  real sparse_iar_lpdf(vector phi, real tau,
    int[,] W_sparse, vector D_sparse, int n, int W_n) {
      row_vector[n] phit_D; // phi' * D
      row_vector[n] phit_W; // phi' * W

      phit_D = (phi .* D_sparse)';
      phit_W = rep_row_vector(0, n);
      for (i in 1:W_n) {
        phit_W[W_sparse[i, 1]] = phit_W[W_sparse[i, 1]] + phi[W_sparse[i, 2]];
        phit_W[W_sparse[i, 2]] = phit_W[W_sparse[i, 2]] + phi[W_sparse[i, 1]];
      }

      return 0.5 * ((n - 1) * log(tau)
                    - tau * (phit_D * phi - (phit_W * phi)));
  }
}

data {
  int<lower = 1> n;                       // sample size
  int<lower = 1> p;                       // number of columns in design matrix
  matrix[n, p] X;                         // design matrix
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
}

parameters {
  real alpha0;
  vector[J] phi1;
  vector[J] phi_diff[n_year - 1];
  real<lower = 0> tau_phi;
  real<lower = 0> tau_phi1;
  real<lower = 0, upper = 1> gamma;
}

transformed parameters {
  vector[J] phi[n_year];
  vector[n] log_lambda;

  phi[1] = phi1 + alpha0;
  for (i in 2:n_year) {
    phi[i] = phi[i - 1] * gamma + phi_diff[i - 1];
  }

  log_lambda = log_offset;
  for (i in 1:n) {
    log_lambda[i] = log_lambda[i] + phi[year[i]][reg[i]];
  }
}

model {
  phi1 ~ sparse_iar(tau_phi1, W_sparse, D_sparse, J, W_n);
  for (i in 1:(n_year-1)) {
    phi_diff[i] ~ sparse_iar(tau_phi, W_sparse, D_sparse, J, W_n);
  }
  alpha0 ~ normal(0, 1);
  tau_phi1 ~ gamma(2, 2);
  tau_phi ~ gamma(2, 2);
  y ~ poisson_log(log_lambda);
}
