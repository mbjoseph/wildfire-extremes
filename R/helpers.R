gev_nll <- function(par, z, tol = 1E-6) {
  # args:
  #   par: parameter vector
  #   z: vector of block maxima
  #   tol: tolerance for when xi is functionally zero
  # returns:
  #   negative log likelihood of the GEV distribution
  mu <- par['mu']
  sigma <- exp(par['logsigma'])
  xi <- par['xi']
  compute_nll(z, mu, sigma, xi, tol)
}


compute_nll <- function(z, mu, sigma, xi, tol) {
  m <- length(z)
  zstar <- (z - mu) / sigma
  xz <- xi * zstar

  if (any(1 + xi * zstar <= 0)) {
    nll <- 1E6
  } else {
    if (abs(xi) < tol) {
      # xi is basically zero, so use Gumbel limit of GEV
      nll <- m * log(sigma) + sum(zstar) + sum(exp(-zstar))
    } else {
      # use GEV
      nll <- m * log(sigma) +
        (1 + 1 / xi) * sum(log(1 + xz)) +
        sum((1 + xz) ^ (-1/xi))
    }
  }
  nll
}
