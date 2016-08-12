source('R/eda.R')
source('R/helpers.R')

# estimate parameters for GEV with timestep-specific maxima
inits <- c(mu = mean(extremes$max_y),
           logsigma = log(sd(extremes$max_y)),
           xi = 0)

estimates <- optim(par = inits,
                   fn = gev_nll,
                   z = extremes$max_y,
                   hessian = TRUE,
                   control = list(maxit = 1E5))
estimates

# compute standard errors
I <- estimates$hessian  # information matrix
Psi <- solve(I)         # covariance matrix
se <- sqrt(diag(Psi))   # assumes the sampling distribution is MvN
se

## using ismev
library(ismev)
mfit <- gev.fit(extremes$max_y)
