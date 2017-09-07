library(scales)
library(gridExtra)
library(spdep)
library(viridis)
library(tidyverse)
library(lubridate)
library(rstan)
library(ggrepel)
library(cowplot)
library(corrplot)
library(sf)
source("R/01-clean_data.R")


# Visualize the distribution of number of fires
count_df %>%
  ggplot(aes(ym, n_fire)) +
  geom_line() +
  facet_wrap(~NA_L3NAME) +
  scale_y_log10()

# Visualize changes in maxima
mtbs %>%
  tbl_df %>%
  group_by(ym, NA_L3NAME) %>%
  summarize(max = max(P_ACRES)) %>%
  ggplot(aes(ym, max)) +
  geom_point() +
  facet_wrap(~ NA_L3NAME)




# Create spatial neighbors ------------------------------------------------
nb <- as(ecoregions, "Spatial") %>%
  poly2nb
W <- nb %>%
  nb2mat(zero.policy = TRUE, style = 'B')
plot(raster(W), col = viridis(10))

# aggregate neighbor matrix to ecoregion level
ecoregion_df <- as(ecoregions, "Spatial") %>%
  data.frame

W_ag <- spdep::aggregate.nb(nb, IDs = ecoregion_df$NA_L3NAME) %>%
  nb2mat(zero.policy = TRUE, style = "B")
plot(raster(W_ag), col = viridis(10))


# visualize spatiotemporal covariates
ecoregion_summaries %>%
  gather(var, val, -NA_L3NAME, -year, -month, -ym) %>%
  ggplot(aes(ym, val)) +
  geom_line() +
  facet_wrap(~ var, scales = "free_y") +
  scale_y_log10()





er_df <- dplyr::select(data.frame(ecoregions), NA_L3NAME, NA_L2NAME, NA_L1NAME) %>%
  as_tibble

st_covs <- ecoregion_summaries %>%
  mutate(dim = 1) %>%
  full_join(mutate(ecoregion_summaries, dim = 2)) %>%
  mutate(cpet = c(scale(pet)),
         cpr = c(scale(pr)),
         ctmx = c(scale(tmmx)),
         cvs = c(scale(vs)),
         dim = factor(dim),
         timestep_factor = factor(ym)) %>%
  left_join(er_df) %>%
  filter(!NA_L2NAME == "UPPER GILA MOUNTAINS (?)")

st_covs <- st_covs[!duplicated(st_covs), ]
st_covs$id <- 1:nrow(st_covs)

# Create design matrices for each timestep
master_idx_df <- filter(st_covs, as.numeric(timestep_factor) == 0) # dummy (0 rows)
n_timestep <- length(levels(st_covs$timestep_factor))

X <- array(dim = c(n_timestep, nrow(st_covs) / n_timestep, 8)) # last dimension is number of cols in design matrix
for (i in 1:n_timestep) {
  master_idx_df <- bind_rows(master_idx_df,
                             filter(st_covs, as.numeric(timestep_factor) == i))
  X[i, , ] <-  filter(st_covs, as.numeric(timestep_factor) == i) %>%
    model.matrix(~ cpr * dim + ctmx * dim + cvs * dim,
                 data = .)
}

# create block diagonal adjacency matrix
A_t <- bdiag(replicate(2, list(W_ag)))
image(A_t)

# compute multivariate basis vectors
r <- 10 # number of basis vectors

S_X <- array(dim = c(n_timestep, nrow(st_covs) / n_timestep, r))
for (i in 1:n_timestep) {
  G <- (diag(nrow(X[1, , ])) -
          X[i, , ] %*% solve(t(X[i, , ]) %*% X[i, , ]) %*% t(X[i, , ])) %*%
    A_t %*%
    (diag(nrow(X[1, , ])) -
       X[i, , ] %*% solve(t(X[i, , ]) %*% X[i, , ]) %*% t(X[i, , ]))
  eG <- eigen(G)
  basis_vectors <- eG$vectors
  S_X[i, , ] <- basis_vectors[, 1:r]
}

# visualize multivariate basis functions
G_unrolled <- S_X[1, , ]
for (i in 2:n_timestep) G_unrolled = rbind(G_unrolled, S_X[i, , ])

# g_df <- tbl_df(G_unrolled) %>%
#   mutate(id = 1:n()) %>%
#   full_join(st_covs) %>%
#   distinct(NA_L3NAME, year, month, dim, .keep_all = TRUE) %>%
#   dplyr::select(NA_L3NAME, year, month, dim, starts_with("V"), -vs) %>%
#   filter(year == 1983, month == 1) %>%
#   dplyr::select(-year, -month) %>%
#   gather(var, value, -NA_L3NAME, -dim)
#
# ecoregions %>%
#   full_join(g_df) %>%
#   ggplot(aes(fill = value)) +
#   geom_sf() +
#   facet_grid(dim ~ var) +
#   scale_fill_gradient2()
#
# ggsave(filename = "fig/st-basis.pdf", width = 13, height = 8)

# confirm that the basis functions are independent of the covariates
tbl_df(G_unrolled) %>%
  mutate(id = 1:n()) %>%
  full_join(st_covs) %>%
  dplyr::select(starts_with("V",
                            ignore.case = FALSE),
                starts_with("c"),
                -cpet) %>%
  cor %>%
  corrplot

# construct the propagator matrix
Bt <- array(dim = c(r, dim(X)[3] + r, n_timestep))
Mt <- array(dim = c(n_timestep, r, r))
for (i in 1:n_timestep) {
  Bt[, , i] <- cbind(t(S_X[i, , ]) %*% X[i, , ], diag(r))
  G_B <- (diag(r) -
            Bt[, , i] %*% MASS::ginv(t(Bt[, , i]) %*% Bt[, , i]) %*% t(Bt[, , i])) %*%
    diag(r) %*%
    (diag(r) -
       Bt[, , i] %*% MASS::ginv(t(Bt[, , i]) %*% Bt[, , i]) %*% t(Bt[, , i]))
  eGB <- eigen(G_B)
  Mt[i, , ] <- eGB$vectors[, 1:r]
}




# simulate basis function coefficients
eta <- matrix(nrow = n_timestep, ncol = r)
R_eta <- rcorrmatrix(r, alphad = 1)
LR_eta <- t(chol(R_eta))
sigma_eta <- .2
sigma_0 <- .5
eta[1, ] <- sigma_0 * c(LR_eta %*% rnorm(r))
for (i in 2:n_timestep) {
  eta[i, ] <- Mt[i, , ] %*% eta[i - 1, ] + sigma_eta * c(LR_eta %*% rnorm(r))
}

# simulate coef for fixed effects
(beta <- c(-5, 3, rnorm(dim(X)[3] - 2, mean = 0, sd = .01)))

# process model
sigma_y <- .1
Y <- matrix(nrow = dim(X)[2], ncol = n_timestep)
for (i in 1:n_timestep) {
  Y[, i] <- S_X[i, , ] %*% eta[i, ] + X[i, , ] %*% beta +
    rnorm(length(Y[, i]), sd = sigma_y)
}

y_df <- Y %>%
  tbl_df %>%
  gather(timestep, z) %>%
  mutate(z = 0.2 * z,
         dim = rep(rep(1:2, each = length(unique(st_covs$NA_L3NAME))), n_timestep),
         dim = factor(dim),
         id = 1:n(),
         timestep = parse_number(timestep),
         timestep_factor = levels(st_covs$timestep_factor)[timestep]) %>%
  bind_cols(dplyr::select(master_idx_df, - dim))



# visualize latent process
simple_ecoregions <- ecoregions %>%
  as('Spatial') %>%
  ms_simplify(keep = 0.01) %>%
  as('sf')


simple_ecoregions %>%
  full_join(y_df) %>%
  filter(timestep < 3) %>%
  ggplot(aes(fill = z)) +
  geom_sf() +
  facet_wrap(dim ~ timestep_factor) +
  scale_fill_gradient2()
ggsave(filename = "fig/us-process-sim.pdf", width = 16, height = 14)


y_df %>%
  ggplot(aes(timestep, z, group = NA_L3NAME)) +
  geom_line(alpha = .1) +
  theme_minimal() +
  facet_wrap(~dim)


# Simulate number of fires in each pixelXyear
n_fire_sims <- y_df %>%
  filter(dim == 1) %>%
  mutate(n_fire = rpois(n(), lambda = exp(z)))

n_fire_sims %>%
  ggplot(aes(exp(z), n_fire)) +
  geom_point()

n_fire_sims %>%
  ggplot(aes(n_fire)) +
  geom_histogram() +
  scale_x_log10()

n_fire_sims %>%
  ggplot(aes(exp(z))) +
  geom_histogram()

# For each fire event, simulate fire size
n_fires <- sum(n_fire_sims$n_fire)



fire_size_df <- tibble(NA_L3NAME = rep(n_fire_sims$NA_L3NAME,
                                       n_fire_sims$n_fire),
                       timestep_factor = rep(n_fire_sims$timestep_factor,
                                  n_fire_sims$n_fire))
fire_size_df <- y_df %>%
  filter(dim == 2) %>%
  dplyr::select(NA_L3NAME, timestep_factor, z) %>%
  mutate(timestep_factor = factor(timestep_factor,
                            levels = levels(master_idx_df$timestep_factor)
  ),
  timestep = as.numeric(timestep_factor)) %>%
  right_join(fire_size_df)

stopifnot(nrow(fire_size_df) == n_fires)

# simulate log fire sizes
fire_size_df <- fire_size_df %>%
  mutate(z_obs = rnorm(n(), mean = z, sd = .1))

fire_size_df %>%
  ggplot(aes(z, z_obs)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")


# Now find indices for each count and each fire size to match to mean matrix
count_idx <- rep(NA, nrow(n_fire_sims))
for (i in seq_along(count_idx)) {
  count_idx[i] <- which(
    master_idx_df$NA_L3NAME == n_fire_sims$NA_L3NAME[i] &
    master_idx_df$timestep_factor == n_fire_sims$timestep_factor[i] &
    master_idx_df$dim == "1")
}

size_idx <- rep(NA, nrow(fire_size_df))
for (i in seq_along(size_idx)) {
  size_idx[i] <-  which(master_idx_df$NA_L3NAME == fire_size_df$NA_L3NAME[i] &
                          master_idx_df$timestep_factor == fire_size_df$timestep_factor[i] &
                          master_idx_df$dim == "2")
}


stan_d <- list(n = length(unique(master_idx_df$NA_L3NAME)),
               T = n_timestep,
               L = 2,
               p = ncol(X[1, ,]),
               r = r,
               X = X,
               S = S_X,
               M = Mt,
               counts = n_fire_sims$n_fire,
               count_idx = count_idx,
               n_fire = n_fires,
               sizes = fire_size_df$z_obs,
               size_idx = size_idx)


# fit model
m_init <- stan_model("stan/mstm.stan", auto_write = TRUE)

m_fit <- sampling(m_init, data = stan_d, cores = 3, chains = 3, iter = 500)
traceplot(m_fit, inc_warmup = TRUE)

traceplot(m_fit, pars = "sigma_size")
traceplot(m_fit, pars = "sigma_eta")

post <- rstan::extract(m_fit)

# Evaluate recovery of spatiotemporal basis coefficients
eta_df <- eta %>%
  reshape2::melt(varnames = c("timestep_factor", "basis_dim"),
                 value.name = "true_eta")


post$eta %>%
  reshape2::melt(varnames = c("iter", "basis_dim", "timestep_factor")) %>%
  tbl_df %>%
  group_by(basis_dim, timestep_factor) %>%
  summarize(med = median(value),
            lo = quantile(value, .025),
            hi = quantile(value, .975)) %>%
  ungroup() %>%
  full_join(eta_df) %>%
  ggplot(aes(true_eta, med)) +
  geom_point() +
  facet_wrap(~basis_dim, scales = "free") +
  geom_segment(aes(xend = true_eta, y = lo, yend = hi)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  xlab("True basis coefficient") +
  ylab("Estimated basis coefficient") +
  ggtitle("Basis coefficient recovery")
ggsave("basis-coef-recovery.pdf", width = 8, height = 6)


# Evaluate recovery of fixed effects
beta_df <- tibble(idx = 1:length(beta),
                  true_value = beta)

post$beta %>%
  reshape2::melt(varnames = c("iter", "idx")) %>%
  full_join(beta_df) %>%
  tbl_df %>%
  ggplot(aes(x = value)) +
  geom_density() +
  facet_wrap(~ paste("Coefficient", idx), scales = "free") +
  geom_vline(aes(xintercept = true_value), linetype = "dashed") +
  xlab("Value") +
  ylab("Posterior density") +
  ggtitle("Recovery of fixed effect coefficients")


post$mu %>%
  reshape2::melt(varnames = c("iter", "idx")) %>%
  mutate(NA_L3NAME = master_idx_df$NA_L3NAME[idx],
         timestep_factor = master_idx_df$timestep_factor[idx],
         dim = master_idx_df$dim[idx]) %>%
  full_join(y_df) %>%
  tbl_df %>%
  group_by(NA_L3NAME, timestep_factor, dim) %>%
  summarize(med = median(value),
            lo = quantile(value, .025),
            hi = quantile(value, .975),
            true = unique(z)) %>%
  ggplot(aes(true, med)) +
  geom_segment(aes(xend = true, y = lo, yend = hi)) +
  facet_grid(dim ~ timestep_factor) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")


