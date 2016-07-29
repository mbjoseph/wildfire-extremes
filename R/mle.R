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

par(mfrow = c(1, 3))
profile_fit(estimates,
            param = 'mu',
            range = c(400, 1000),
            z = extremes$max_y)

profile_fit(estimates,
            param = 'logsigma',
            range = c(5.5, 7.5),
            z = extremes$max_y)

legend('top',
       lty = c(1, 2, 0),
       pch = c(NA, NA, 1),
       col = c('black', 'blue', 'red'),
       legend = c('NLL profile', '95% CI threshold', 'MLE'),
       bty = 'n')


profile_fit(estimates,
            param = 'xi',
            range = c(-.3, .7),
            z = extremes$max_y)




## using ismev
library(ismev)
mfit <- gev.fit(extremes$max_y)
