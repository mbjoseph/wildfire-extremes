# simulating block maxima from lognormal distribution with varying sample size
par(mfrow = c(1, 1))

# view large sample from lognormal distribution
hist(exp(rnorm(1E5)), breaks = 100, freq = FALSE)

# simulate sample maxima with varying sample size
n <- 1:1E3
ymax <- unlist(lapply(X = n, FUN = function(x) max(exp(rnorm(x, 1, 3)))))

# what is the relationship between sample size and the sample max?
plot(n, ymax)
plot(n, log(ymax))
plot(log(n), log(ymax))
