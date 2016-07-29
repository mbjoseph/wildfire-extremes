/**
 * Log probability density of generalized extreme value (GEV) distribution
 *
 * @param y Observation of an extreme value (e.g., a block maximum)
 * @param mu Location parameter
 * @param sigma Scale parameter
 * @param xi Shape parameter
*/

functions {
  real gev_lpdf(real y, real mu, real sigma, real xi) {
    real ystar;
    real xy_p1;

    ystar = (y - mu) / sigma;
    xy_p1 = xi * ystar + 1;
    return -log(sigma) - (1 + 1 / xi) * log(xy_p1) - (xy_p1) ^ (-1 / xi);
  }
}

data {
  int n;
  vector[n] y;
}

parameters {
  real mu;
  real<lower = 0> sigma;
  real xi;
}

model {
  xi ~ normal(0, 10);

  for (i in 1:n) {
    y[i] ~ gev(mu, sigma, xi);
  }
}
