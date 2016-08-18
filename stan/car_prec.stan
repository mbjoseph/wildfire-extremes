data {
  int n;
  int p;
  int y[n];
  vector[n] log_offset;
  matrix[n, p] X;
  matrix[n, n] D;
  matrix[n, n] W;
}
transformed data {
  vector[n] zeros;
  zeros = rep_vector(0, n);
}
parameters {
  real<lower = 0> tau;
  real<lower = 0, upper = 1> alpha;
  vector[p] beta;
  vector[n] phi;
}
model {
  tau ~ normal(0, 5);
  beta ~ normal(0, 10);
  phi ~ multi_normal_prec(zeros, tau * (D - alpha * W));
  y ~ poisson_log(X * beta + phi + log_offset);
}
