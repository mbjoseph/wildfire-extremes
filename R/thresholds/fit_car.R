
# CAR models --------------------------------------------------------------
source('R/explore.R')

# simulate realizations from CAR prior
D <- diag(rowSums(W))
invD <- diag(1 / diag(D))
B <- invD %*% W
plot(raster(B))

alpha <- .5
n <- nrow(W)
sigma <- 1.5
sig_sq <- sigma ** 2
tau <- 1 / sig_sq
Tau <- (tau * D) %*% (diag(n) - alpha * B)
Sigma <- solve(Tau)
plot(raster(Sigma))
L <- t(chol(Sigma))
phi <- L %*% rnorm(n)

sim_df <- data.frame(us_l3name = rownames(W),
                     phi = phi)
sim_df$area <- er_df %>%
  dplyr::select(us_l3name, area) %>%
  group_by(us_l3name) %>%
  summarize(area = unique(area)) %>%
  dplyr::select(area) %>%
  unlist()


sim_df$y <- rpois(n, exp(-23 + sim_df$phi + log(sim_df$area)))
hist(sim_df$y)
car_df <- full_join(er_df, sim_df)

ggplot(car_df, aes(long, lat, group = group)) +
  geom_polygon(aes(fill = y), color = NA) +
  coord_equal() +
  labs(x = "Longitude", y = "Latitude") +
  scale_fill_gradient(low = 'black', high = 'orange') +
  theme_map

agg_d <- car_df %>%
  group_by(us_l3name) %>%
  summarize(area = unique(area), n_fires = unique(n_fires), y = unique(y))

library(rstan)
options(mc.cores = parallel::detectCores())
stan_d <- list(n = n,
               X = model.matrix(~1, data = agg_d),
               p = 1,
               D = D,
               W = W,
               y = agg_d$y,
               log_offset = log(agg_d$area))
m_init <- stan('stan/car_prec.stan', data = stan_d, iter = 1, chains = 1)
m_fit <- stan(fit = m_init, data = stan_d, iter = 4000, chains = 4)
m_fit
traceplot(m_fit)

# check parameter recovery
post <- extract(m_fit)
phi_est <- apply(post$phi, 2, function(x) {
  c(lo = quantile(x, .1),
    med = median(x),
    hi = quantile(x, .9))
}) %>%
  t() %>%
  as.data.frame()
phi_est$true_phi <- c(phi)

ggplot(phi_est, aes(x = true_phi, y = med)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
  geom_segment(aes(y = `lo.10%`, yend = `hi.90%`, xend = true_phi))


# Sparse CAR --------------------------------------------------------------
D_sparse <- diag(D)
# W is two vectors corresponding to the non-zero pairs
W_sparse <- which(W == 1, arr.ind = TRUE)
W_sparse <- W_sparse[W_sparse[,1] < W_sparse[,2],]  # remove duplicates (because matrix is symmetric)
W_n <- dim(W_sparse)[1]
W1 <- W_sparse[,1]
W2 <- W_sparse[,2]

# so, about that determinant...
# get eigenvalues for D^(-.5) %*% W %*% D^(-.5)
isqrtD <- diag(1/(sqrt(diag(D))))
lambda <- eigen(isqrtD %*% W %*% isqrtD)$values

stan_d <- list(n = n,
               X = model.matrix(~1, data = agg_d),
               p = 1,
               D = D,
               W = W,
               y = agg_d$y,
               log_offset = log(agg_d$area),
               W_n = W_n,
               W1 = W1,
               W2 = W2,
               D_sparse = D_sparse,
               lambda = lambda)

sp_init <- stan('stan/sparse_car.stan', data = stan_d, iter = 1, chains = 1)
sp_fit <- stan(fit = sp_init, data = stan_d, iter = 2000, chains = 4)
sp_fit
traceplot(sp_fit, pars = c('alpha', 'beta', 'tau', 'phi[1]', 'phi[2]'))

post <- rstan::extract(sp_fit)
str(post)
phi_est2 <- apply(post$phi, 2, function(x) {
  c(lo = quantile(x, .1),
    med = median(x),
    hi = quantile(x, .9))
}) %>%
  t() %>%
  as.data.frame()
phi_est2$true_phi <- c(phi)

ggplot(phi_est2, aes(x = true_phi, y = med)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
  geom_segment(aes(y = `lo.10%`, yend = `hi.90%`, xend = true_phi))


# join and compare
phi_est$model <- 'full'
phi_est2$model <- 'sparse'
phi_df <- full_join(phi_est, phi_est2)

ggplot(phi_df, aes(x = true_phi, y = med, color = model)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
  geom_segment(aes(y = `lo.10%`, yend = `hi.90%`, xend = true_phi))


# Fit areal log gaussian cox process model using actual fire freq. data
agg_d$n_fires[is.na(agg_d$n_fires)] <- 0
stan_d <- list(n = n,
               X = model.matrix(~1, data = agg_d),
               p = 1,
               D = D,
               W = W,
               y = agg_d$n_fires,
               log_offset = log(agg_d$area),
               W_n = W_n,
               W1 = W1,
               W2 = W2,
               D_sparse = D_sparse,
               lambda = lambda)
sp_fit <- stan(fit = sp_init, data = stan_d, iter = 5000, chains = 4,
               warmup = 2000, pars = c('log_mu', 'beta', 'phi', 'tau', 'alpha'),
               control = list(max_treedepth = 13))
sp_fit
traceplot(sp_fit, pars = c('alpha', 'tau', 'beta', 'phi[1]', 'phi[2]', 'phi[3]', 'lp__'))
pairs(sp_fit, pars = c('alpha', 'tau', 'beta', 'phi[1]', 'phi[2]', 'phi[3]', 'lp__'))
plot(sp_fit, pars = 'phi')
plot(sp_fit, pars = 'log_mu')

# map the results
post <- rstan::extract(sp_fit)
post$mu <- exp(post$log_mu)
res_df <- data.frame(mu_med = apply(post$mu, 2, median),
                     mu_lo = apply(post$mu, 2, quantile, 0.025),
                     mu_hi = apply(post$mu, 2, quantile, 0.975),
                     mu_se = apply(post$mu, 2, sd),
                     phi_med = apply(post$phi, 2, median),
                     phi_se = apply(post$phi, 2, sd),
                     us_l3name = rownames(W))

full_join(res_df, car_df) %>%
  ggplot(aes(long, lat, group = group)) +
  geom_polygon(aes(fill = mu_med), color = NA) +
  coord_equal() +
  labs(x = "Longitude", y = "Latitude") +
  scale_fill_gradient(low = 'black', high = 'orange') +
  theme_map

full_join(res_df, car_df) %>%
  ggplot(aes(long, lat, group = group)) +
  geom_polygon(aes(fill = mu_se), color = NA) +
  coord_equal() +
  labs(x = "Longitude", y = "Latitude") +
  scale_fill_gradient(low = 'black', high = 'orange') +
  theme_map

full_join(res_df, car_df) %>%
  ggplot(aes(long, lat, group = group)) +
  geom_polygon(aes(fill = phi_med), color = NA) +
  coord_equal() +
  labs(x = "Longitude", y = "Latitude") +
  scale_fill_gradient(low = 'black', high = 'orange') +
  theme_map

full_join(res_df, car_df) %>%
  ggplot(aes(long, lat, group = group)) +
  geom_polygon(aes(fill = phi_se), color = NA) +
  coord_equal() +
  labs(x = "Longitude", y = "Latitude") +
  scale_fill_gradient(low = 'black', high = 'orange') +
  theme_map


## How well do the estimates match the empirical counts?
full_join(res_df, agg_d) %>%
  group_by(us_l3name) %>%
  summarize(mu_med = unique(mu_med),
            mu_lo = unique(mu_lo),
            mu_hi = unique(mu_hi),
            n_fires = n_fires) %>%
  ggplot(aes(x = n_fires)) +
  geom_point(aes(y = mu_med)) +
  geom_segment(aes(xend = n_fires, y = mu_lo, yend = mu_hi)) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed')

