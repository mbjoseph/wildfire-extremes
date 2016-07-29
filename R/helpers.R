gev_nll <- function(par, z, tol = 1E-6) {
  # args:
  #   par: parameter vector (length 3: mu, log(sigma), xi)
  #   y: block maxima
  #   tol: tolerance for when xi is functionally zero
  # returns:
  #   negative log likelihood of the GEV distribution

  mu <- par['mu']
  sigma <- exp(par['logsigma'])
  xi <- par['xi']

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


# generate likelihood profiles
profile_fit <- function(fit, param, range, z, n = 1000) {
  par_range <- seq(range[1], range[2], length.out = n)
  nll_out <- rep(NA, n)
  mle <- fit$par
  for (i in seq_along(par_range)) {
    par <- mle
    par[param] <- par_range[i]
    nll_out[i] <- gev_nll(par, z)
  }
  cbind(par_range, nll_out) %>%
    plot(xlab = param, ylab = 'Negative log likelihood', type = 'l')
  points(x = mle[param], y = fit$value, col = 'red')
  # draw line to delineate likelihood ratio based confidence intervals
  abline(h = fit$value + qchisq(0.95, df = 1) / 2, lty = 2, col = 'blue')
}
