data {
  int<lower = 1> N; // # spatial units
  int<lower = 1> T; // # timesteps
  int<lower = 1> L; // # dimensions
  int<lower = 1> p; // # columns in design matrix
  int<lower = 1> p_c; // # cols in count X
  vector[N * T] log_area;
  int<lower = 0> counts[N * T]; // # of events in each spatial unit * timestep

  // sparse matrix for counts
  int<lower = 1> n_wc;
  vector[n_wc] wc;
  int<lower = 1> vc[n_wc];
  int<lower = 1> uc[N * T + 1];

  int n_fire;
  vector[n_fire] sizes;

  // sparse matrix for fire sizes
  int<lower = 1> n_w;
  vector[n_w] w;
  int<lower = 1> v[n_w];
  int<lower = 1> u[n_fire + 1];
}

parameters {
  vector[p] betaR;
  vector<lower = 0>[p] lambda;
  real<lower = 0> tau;

  vector[p_c] betaR_c;
  vector<lower = 0>[p_c] lambda_c;
  real<lower = 0> tau_c;

  vector[2] alpha; // intercept
  real<lower = 0> sigma_size;
  vector[N * T] count_epsR;
  real<lower = 0> sigma_eps;

  vector<lower = 0>[2] c_sq;
}

transformed parameters {
  vector[n_fire] mu;
  vector[N * T] mu_counts;
  vector[p] beta;
  vector[p] lambda_tilde;
  vector[p_c] beta_c;
  vector[p_c] lambda_tilde_c;

  // regularized horseshoe prior
  lambda_tilde = sqrt(c_sq[1] * square(lambda) ./
                      (c_sq[1] + tau^2 * square(lambda)));
  lambda_tilde_c = sqrt(c_sq[2] * square(lambda_c) ./
                        (c_sq[2] + tau_c^2 * square(lambda_c)));

  beta = betaR .* lambda_tilde * tau;
  beta_c = betaR_c .* lambda_tilde_c * tau_c;

  mu = alpha[1] + csr_matrix_times_vector(n_fire, p, w, v, u, beta);

  mu_counts = alpha[2]
                + csr_matrix_times_vector(N*T, p_c, wc, vc, uc, beta_c)
                + count_epsR * sigma_eps
                + log_area;
}

model {
  count_epsR ~ normal(0, 1);
  sigma_eps ~ normal(0, 1);
  sigma_size ~ normal(0, 1);
  alpha ~ normal(0, 1);

  c_sq[1] ~ inv_gamma(2, 2);
  c_sq[2] ~ inv_gamma(1, 1);

  betaR ~ normal(0, 1);
  lambda ~ cauchy(0, 1);
  tau ~ student_t(5, 0, 1);

  betaR_c ~ normal(0, 1);
  lambda_c ~ cauchy(0, 1);
  tau_c ~ student_t(5, 0, 1);

  // number of fires
  counts ~ poisson_log(mu_counts);

  // fire sizes
  sizes ~ normal(mu, sigma_size);
}
