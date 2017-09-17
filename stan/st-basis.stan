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
  int<lower = 1> n_u;
  int<lower = 1> u[n_u]l

  int<lower = 0> counts[N * T]; // # of events in each spatial unit * timestep
  int<lower = 1, upper = nrowX> count_idx[N * T]; // which row in X for each count
  int n_fire;
  vector[n_fire] sizes;
  int<lower = 1, upper = nrowX> size_idx[n_fire];
}

parameters {
  vector[p] beta;
  real<lower = 0> sigma_size;
}

transformed parameters {
  vector[nrowX] mu;
  mu = csr_matrix_times_vector(nrowX, p, w, v, u, beta);
}

model {
  beta ~ normal(0, 5);

  // number of fires
  counts ~ poisson_log(mu[count_idx]);

  // fire sizes
  sigma_size ~ normal(0, 1);
  sizes ~ normal(mu[size_idx], sigma_size);
}
