data {
  int n;
  int p;
  int y[n];
  matrix[n, p] X;
  matrix[n, n] B;
  matrix[n, n] D;
  matrix[n, n] I;
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
  beta ~ normal(0, 1);
  phi ~ multi_normal_prec(zeros, tau * D * (I - alpha * B));
  y ~ poisson_log(X * beta + phi);
}
