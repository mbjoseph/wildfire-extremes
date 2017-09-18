data {
  int<lower = 1> N; // # spatial units
  int<lower = 1> T; // # timesteps
  int<lower = 1> L; // # dimensions
  int<lower = 1> p; // # columns in design matrix
  int<lower = 1, upper = N * L * T> nrowX; // # rows in design matrix to use

  // sparse representation of X
  int<lower = 1> n_w;
  vector[n_w] w;
  int<lower = 1> v[n_w];
  int<lower = 1> u[nrowX + 1];

  int<lower = 0> counts[N * T]; // # of events in each spatial unit * timestep
  int<lower = 1, upper = nrowX> count_idx[N * T]; // which row in X for each count
  int n_fire;
  vector[n_fire] sizes;
  int<lower = 1, upper = nrowX> size_idx[n_fire];
  vector[N * T] log_area;
}

parameters {
  vector[p] betaR;
  vector<lower = 0>[p] lambda;
  real<lower = 0> tau;
  real<lower = 0> sigma_size;
}

transformed parameters {
  vector[nrowX] mu;
  vector[N * T] mu_counts;
  vector[p] beta;

  beta = betaR .* lambda * tau;
  mu = csr_matrix_times_vector(nrowX, p, w, v, u, beta);
  mu_counts = mu[count_idx] + log_area;
}

model {
  betaR ~ normal(0, 1);
  lambda ~ cauchy(0, 1);
  tau ~ normal(0, 1);

  // number of fires
  counts ~ poisson_log(mu_counts);

  // fire sizes
  sigma_size ~ normal(0, 1);
  sizes ~ normal(mu[size_idx], sigma_size);
}
