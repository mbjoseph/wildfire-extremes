functions {
  /**
  * Return the log probability of a proper scale free conditional autoregressive
  * (CAR) prior with a sparse representation for the adjacency matrix
  *
  * @param phi Vector containing the parameters with a CAR prior
  * @param alpha Dependence (usually spatial) parameter for the CAR prior (real)
  * @param W_sparse Sparse representation of adjacency matrix (int array)
  * @param n Length of phi (int)
  * @param W_n Number of adjacent pairs (int)
  * @param D_sparse Number of neighbors for each location (vector)
  * @param lambda Eigenvalues of D^{-1/2}*W*D^{-1/2} (vector)
  *
  * @return Log probability density of CAR prior up to additive constant
  */
  real sparse_car_lpdf(vector phi, real alpha,
    int[,] W_sparse, vector D_sparse, vector lambda, int n, int W_n) {
      row_vector[n] phit_D; // phi' * D
      row_vector[n] phit_W; // phi' * W
      vector[n] ldet_terms;

      phit_D = (phi .* D_sparse)';
      phit_W = rep_row_vector(0, n);
      for (i in 1:W_n) {
        phit_W[W_sparse[i, 1]] = phit_W[W_sparse[i, 1]] + phi[W_sparse[i, 2]];
        phit_W[W_sparse[i, 2]] = phit_W[W_sparse[i, 2]] + phi[W_sparse[i, 1]];
      }

      for (i in 1:n) ldet_terms[i] = log1m(alpha * lambda[i]);
      return 0.5 * (sum(ldet_terms)
                    - (phit_D * phi - alpha * (phit_W * phi)));
  }
}
data {
  // fire counts
  int<lower = 1> n;                       // sample size
  int<lower = 0> y[n];                    // number of fires
  int<lower = 1> J;                       // number of ecoregions
  int<lower = 1, upper = J> reg[n];       // ecoregion index
  vector[J] log_offset;                   // log area
  matrix<lower = 0, upper = 1>[J, J] W;   // adjacency matrix
  int W_n;                                // number of adjacent region pairs
  int n_year;
  int<lower = 1, upper = n_year> year[n];

  // fire sizes
  int<lower = 1> nz;  // number of fire size observations
  vector[nz] z;       // log fire sizes
  int<lower = 1, upper = J> reg_z[nz];       // ecoregion index
  int<lower = 1, upper = n_year> year_z[nz];

  // precip
  vector[nz] prcp_z;
  vector[n] prcp_y;
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
  vector<lower = 0>[2] sd_sd;
  vector[J] gammaR[2];
  vector<lower =0>[2] sd_gamma;
  vector[2] mu_gamma;
  vector<lower = 0, upper = 1>[7] alpha;

  vector[J] phi_yR[n_year];
  vector<lower = 0>[n_year] sigma_phi_y;

  real z0;

  // regional variance
  vector<lower = 0>[J] sigma_z;
  real mu_sigma_z;
  real<lower = 0> sigma_sigma_z;

  vector[J] phi_zR[n_year];
  vector<lower = 0>[n_year] sigma_phi_z;
  vector[J] beta_phi_yR;
  real mu_beta_phi;
  real<lower = 0> sigma_beta_phi;

  vector[J] beta_prcpR[2];
  vector<lower = 0>[2] sigma_prcp;
  vector[2] mu_prcp;
}

transformed parameters {
  matrix[J, n_year] phi_y;
  matrix[J, n_year] phi_z;
  vector[n] log_lambda;
  vector[nz] mu_z;
  vector[J] beta_phi_y;
  vector[J] gamma[2];
  vector[J] beta_prcp[2];

  for (i in 1:2) {
    gamma[i] = gammaR[i] * sd_gamma[i] + mu_gamma[i];
    beta_prcp[i] = beta_prcpR[i] * sigma_prcp[i] + mu_prcp[i];
  }

  beta_phi_y = beta_phi_yR * sigma_beta_phi + mu_beta_phi;

  phi_y[, 1] = phi_yR[1] * sigma_phi_y[1];
  phi_z[, 1] = phi_zR[1] * sigma_phi_z[1];
  for (i in 2:n_year) {
    phi_y[, i] = gamma[1] .* phi_y[, i - 1] + sigma_phi_y[i] * phi_yR[i];
    phi_z[, i] = gamma[2] .* phi_z[, i - 1] + sigma_phi_z[i] * phi_zR[i]
                  + beta_phi_y .* phi_y[, i];
  }

  for (i in 1:n) {
    log_lambda[i] = log_offset[reg[i]] + phi_y[reg[i], year[i]];
  }
  log_lambda = log_lambda + beta_prcp[1][reg] .* prcp_y;

  for (i in 1:nz) {
    mu_z[i] = phi_z[reg_z[i], year_z[i]];
  }
  mu_z = mu_z + z0 + beta_prcp[2][reg_z] .* prcp_z;
}

model {
  for (i in 1:n_year) {
    phi_yR[i] ~ sparse_car(alpha[1], W_sparse, D_sparse, lambda, J, W_n);
    phi_zR[i] ~ sparse_car(alpha[2], W_sparse, D_sparse, lambda, J, W_n);
  }

  sd_sd ~ normal(0, 1);
  sigma_phi_y ~ normal(0, sd_sd[1]);
  sigma_phi_z ~ normal(0, sd_sd[2]);

  z0 ~ normal(0, 5);
  beta_phi_yR ~ sparse_car(alpha[3], W_sparse, D_sparse, lambda, J, W_n);
  mu_beta_phi ~ normal(0, 1);
  sigma_beta_phi ~ lognormal(0, 1);
  sigma_z ~ lognormal(mu_sigma_z, sigma_sigma_z);
  mu_sigma_z ~ normal(0, 1);
  sigma_sigma_z ~ normal(0, .1);

  for (i in 1:2) {
    gammaR[i] ~ sparse_car(alpha[3 + i], W_sparse, D_sparse, lambda, J, W_n);
  }
  mu_gamma ~ normal(0, 1);
  sd_gamma ~ normal(0, 1);

  for (i in 1:2) beta_prcpR[i] ~ sparse_car(alpha[5 + i], W_sparse, D_sparse, lambda, J, W_n);
  sigma_prcp ~ normal(0, 1);
  mu_prcp ~ normal(0, 1);




  y ~ poisson_log(log_lambda);

  z ~ normal(mu_z, sigma_z[reg_z]);
}
//
// generated quantities {
//   matrix[J, n_year] log_lambda_new;
//   matrix[J, n_year] mu_z_new;
//   int y_new[J, n_year];
//   matrix[J, n_year] zmax_new;
//   matrix[J, n_year] area_burned;
//
//   // draw parameters for each ecoregion X year
//   for (i in 1:n_year) {
//     for (j in 1:J) {
//       area_burned[j, i] = 0;
//       zmax_new[j, i] = 0;
//       log_lambda_new[j, i] = log_offset[j] + phi_y[j, i];
//       mu_z_new[j, i] = phi_z[j, i] + z0;
//       y_new[j, i] = poisson_log_rng(log_lambda_new[j, i]);
//       if (y_new[j, i] > 0) {
//       { // simulate fire sizes
//         vector[y_new[j, i]] z_new;
//         for (k in 1:y_new[j, i]) z_new[k] = normal_rng(mu_z_new[j, i], sigma_z[j]);
//         zmax_new[j, i] = max(z_new);
//         area_burned[j, i] = sum(z_new);
//       }
//       }
//     }
//   }
// }
