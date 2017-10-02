data {
  int<lower = 1> N; // # spatial units
  int<lower = 1> T; // # timesteps
  int<lower = 1> p; // # columns in design matrix
  vector[N * T] log_area;

  int<lower = 1> n_count; // number of counts in training set
  int<lower = 0> counts[n_count]; // # of events in each spatial unit * timestep
  int<lower = 1, upper = N * T> count_idx[n_count];

  int n_fire;
  vector[n_fire] sizes;
  int<lower = 1, upper = N * T> burn_idx[n_fire];

  // sparse matrix for fire sizes
  int<lower = 1> n_w;
  vector[n_w] w;
  int<lower = 1> v[n_w];
  int<lower = 1> u[N * T + 1];
}

parameters {
  vector[p] betaR[3];
  vector<lower = 0>[p] lambda[3];
  vector<lower = 0>[3] tau;
  vector[3] alpha; // intercept
  vector[N * T] epsR[3];
  vector<lower = 0>[3] sigma;
  vector<lower = 0>[3] c_sq;
}

transformed parameters {
  vector[N * T] mu[3];
  vector[p] beta[3];
  vector[p] lambda_tilde[3];
  vector[p] lambda_sq[3];

  // regularized horseshoe prior
  for (i in 1:3) {
    lambda_sq[i] = square(lambda[i]);
    lambda_tilde[i] = sqrt(
                            c_sq[i] * lambda_sq[i] ./
                            (c_sq[i] + tau[i]^2 * lambda_sq[i])
                          );
    beta[i] = betaR[i] .* lambda_tilde[i] * tau[i];
    mu[i] = alpha[i] +
              csr_matrix_times_vector(N * T, p, w, v, u, beta[i]) +
              epsR[i] * sigma[i];
  }

  // expected log counts need an offset for area
  mu[3] = mu[3] + log_area;
}

model {
  for (i in 1:3) {
    betaR[i] ~ normal(0, 1);
    lambda[i] ~ cauchy(0, 1);
    epsR[i] ~ normal(0, 1);
  }

  c_sq ~ inv_gamma(1, 1);
  sigma ~ normal(0, 1);
  alpha ~ normal(0, 1);
  tau ~ normal(0, 1);

  // number of fires
  counts ~ poisson_log(mu[3][count_idx]);

  // fire sizes
  sizes ~ weibull(exp(mu[1])[burn_idx], exp(mu[2])[burn_idx]);
}
