
# Verifying variance of ZI Poisson ----------------------------------------

n_iter <- 1000
n <- 10000
p <- .5
lambda <- 4

mu_nb <- 4
phi_nb <- 2


mu <- lambda * p
V <- mu * p + p * (1- p) * (lambda^2 + lambda)

mu_zinb <- mu_nb * p
V_zinb <- p^2 * (mu_nb + mu_nb^2 / phi_nb) +
  mu_nb^2 * p * (1-p) +
  (mu_nb + mu_nb^2 / phi_nb) * p * (1 - p)

means <- rep(NA, n_iter)
vars <- rep(NA, n_iter)

means_zinb <- rep(NA, n_iter)
vars_zinb <- rep(NA, n_iter)


for (i in 1:n_iter) {
  z <- rbinom(n, size = 1, prob = p)
  y <- rpois(n, lambda = lambda)
  y_nb <- rnbinom(n, size = phi_nb, mu = mu_nb)
  x <- z * y
  x_zinb <- z * y_nb
  means[i] <- mean(x)
  vars[i] <- var(x)
  means_zinb[i] <- mean(x_zinb)
  vars_zinb[i] <- var(x_zinb)
}

par(mfrow=c(2, 2))
plot(sort(means))
abline(h = mu, col = 2)
plot(sort(vars))
abline(h = V, col = 2)

plot(sort(means_zinb))
abline(h = mu_zinb, col = 2)
plot(sort(vars_zinb))
abline(h = V_zinb, col = 2)
