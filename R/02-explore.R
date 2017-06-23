library(scales)
library(gridExtra)
library(spdep)
library(viridis)
library(tidyverse)
library(lubridate)
library(rstan)
library(ggrepel)
library(cowplot)

source("R/01-clean_data.R")


# subset for sims
d <- d %>%
  filter(fire_year < 1995)


# Visualize yearly fire size distributions ----------------------------
yd <- d %>%
  group_by(fire_year) %>%
  summarize(mean_size = mean(fire_size),
            max_size = max(fire_size),
            n = n()) %>%
  ungroup %>%
  tidyr::complete(fire_year,
           fill = list(mean_size = NA,
                       max_size = NA)) %>%
  gather(Measure, value, -fire_year)


yd %>%
  ggplot(aes(x = fire_year, y = value)) +
  geom_path() +
  facet_wrap(~ Measure, scales = "free_y") +
  xlab("Time") +
  ylab("Value") +
  theme_minimal()
ggsave("fig/explore/mean-max-ts.pdf")


p1 <- yd %>%
  filter(Measure != "mean_size") %>%
  mutate(Measure = ifelse(Measure == "max_size", "Maximum burn area",
                          "Number of fires > 1000 acres")) %>%
  ggplot(aes(x = fire_year, y = value)) +
  geom_path() +
  facet_wrap(~ Measure, scales = "free_y") +
  xlab("Time") +
  ylab("Value") +
  theme_minimal()

p2 <- yd %>%
  filter(Measure != "mean_size") %>%
  mutate(Measure = ifelse(Measure == "max_size", "Maximum burn area",
                          "Number of fires > 1000 acres")) %>%
  spread(Measure, value) %>%
  ggplot(aes(`Number of fires > 1000 acres`, `Maximum burn area`)) +
  geom_point() +
  theme_minimal() +
  geom_text_repel(aes(label = fire_year))

plot_grid(p1, p2, nrow = 2)
ggsave("fig/explore/maxima-n.pdf", width = 9, height = 5)

# Create spatial neighbors ------------------------------------------------
W <- coarse_rp %>%
  rasterToPolygons %>%
  poly2nb %>%
  nb2mat(zero.policy = TRUE, style = 'B')
plot(raster(W))

plot(rasterToPolygons(coarse_rp))
coarse_rp %>%
  rasterToPolygons %>%
  poly2nb %>%
  plot(col = "red",
       coords = coordinates(rasterToPolygons(coarse_rp)))

# Create spatiotemporal covariate data frame
st_cov_chunk <- coarse_rp %>%
  as.data.frame(xy = TRUE) %>%
  tbl_df %>%
  mutate(pixel_idx = 1:n()) %>%
  gather(variable, value, -x, -y, -pixel_idx) %>%
  separate(variable, c("variable", "year"), sep = "_") %>%
  arrange(pixel_idx, year, variable) %>%
  filter(!is.na(value)) %>%
  spread(variable, value) %>%
  arrange(pixel_idx, year)

st_covs <- st_cov_chunk %>%
  mutate(dim = 1) %>%
  full_join(mutate(st_cov_chunk, dim = 2)) %>%
  mutate(chuss = c(scale(huss)),
         cprec = c(scale(pr)),
         ctmax = c(scale(tasmax)),
         ctmin = c(scale(tasmin)),
         cwas = c(scale(was)),
         year = parse_number(year),
         dim = factor(dim))


# Create design matrices for each timestep

st_covs <- st_covs %>%
  mutate(year_idx = factor(year))

master_idx_df <- filter(st_covs, as.numeric(year_idx) == 0)
X <- array(dim = c(n_year, nrow(st_covs) / n_year, 7)) # 2: bivariate
for (i in 1:n_year) {
  master_idx_df <- bind_rows(master_idx_df, filter(st_covs, as.numeric(year_idx) == i))
  X[i, , ] <-  filter(st_covs, as.numeric(year_idx) == i) %>%
    model.matrix(~ dim * chuss + cprec + ctmax + cwas, data = .)
}

A_t <- bdiag(replicate(2, list(W)))
image(A_t)

# compute multivariate basis vectors
r <- 50 # number of basis vectors

S_X <- array(dim = c(n_year, nrow(st_covs) / n_year, r))
for (i in 1:n_year) {
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
for (i in 2:n_year) G_unrolled = rbind(G_unrolled, S_X[i, , ])

tbl_df(G_unrolled) %>%
  mutate(pixel_idx = rep(rep(unique(st_covs$pixel_idx), times = 2), n_year),
         year = rep(unique(st_covs$year), each = nrow(st_covs) / n_year),
         dim = rep(rep(1:2, each = length(unique(st_covs$pixel_idx))), n_year),
         dim = factor(dim)) %>%
  full_join(st_covs) %>%
  dplyr::select(x, y, year, dim, starts_with("V")) %>%
  gather(basis_dim, value, -x, -y, -year, -dim) %>%
  filter(year == 1984) %>%
  mutate(basis_dim = parse_number(basis_dim),
         Dimension = ifelse(dim == 1, "# Fires", "Burn area"),
         Basis_dim = paste("Basis", sprintf("%02s", basis_dim))) %>%
  ggplot(aes(x, y, fill = value)) +
  geom_tile() +
  facet_wrap(Basis_dim ~ Dimension, ncol = 10) +
  scale_fill_gradient2("Basis value") +
  coord_equal() +
  theme_dark() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("Sample of spatiotemporal basis functions")

ggsave(filename = "fig/st-basis.pdf", width = 13, height = 8)

# confirm that the basis functions are independent of the covariates
tbl_df(G_unrolled) %>%
  bind_cols(st_covs) %>%
  dplyr::select(starts_with("V"), chuss, cprec, ctmax, cwas) %>%
  cor %>%
  image

# construct the propagator matrix
Bt <- array(dim = c(r, dim(X)[3] + r, n_year))
Mt <- array(dim = c(n_year, r, r))
for (i in 1:n_year) {
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
eta <- matrix(nrow = n_year, ncol = r)
R_eta <- rcorrmatrix(r, alphad = 1)
LR_eta <- t(chol(R_eta))
sigma_eta <- .1
sigma_0 <- .1
eta[1, ] <- sigma_0 * c(LR_eta %*% rnorm(r))
for (i in 2:n_year) {
  eta[i, ] <- Mt[i, , ] %*% eta[i - 1, ] + sigma_eta * c(LR_eta %*% rnorm(r))
}

# simulate coef for fixed effects
(beta <- c(-.5, 0, rnorm(dim(X)[3] - 2, mean = 0, sd = .1)))

# process model
Y <- matrix(nrow = dim(X)[2], ncol = n_year)
for (i in 1:n_year) {
  Y[, i] <- S_X[i, , ] %*% eta[i, ] + X[i, , ] %*% beta
}

# visualize latent process
Y %>%
  tbl_df %>%
  gather(sim_year, z) %>%
  mutate(pixel_idx = rep(rep(unique(st_covs$pixel_idx), times = 2), n_year),
         year = rep(unique(st_covs$year), each = nrow(st_covs) / n_year),
         dim = rep(rep(1:2, each = length(unique(st_covs$pixel_idx))), n_year),
         dim = factor(dim)) %>%
  full_join(st_covs) %>%
  ggplot(aes(x, y, fill = z)) +
  geom_tile() +
  facet_grid(dim ~ year) +
  scale_fill_viridis() +
  coord_equal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  ggtitle("Simulated spatiotemporal process")
ggsave(filename = "fig/us-process-sim.pdf", width = 24, height = 5)


Y %>%
  tbl_df %>%
  gather(sim_year, z) %>%
  mutate(pixel_idx = rep(rep(unique(st_covs$pixel_idx), times = 2), n_year),
         year = rep(unique(st_covs$year), each = nrow(st_covs) / n_year),
         dim = rep(rep(1:2, each = length(unique(st_covs$pixel_idx))), n_year),
         dim = factor(dim)) %>%
  full_join(st_covs) %>%
  ggplot(aes(year, z, group = pixel_idx)) +
  geom_line(alpha = .1) +
  theme_minimal() +
  facet_wrap(~dim)

# Simulate number of fires in each pixelXyear
n_fire_sims <- Y %>%
  tbl_df %>%
  gather(sim_year, z) %>%
  bind_cols(st_covs) %>%
  filter(dim == 1) %>%
  mutate(n_fire = rpois(n(), lambda = exp(z)))

n_fire_sims %>%
  ggplot(aes(exp(z), n_fire)) +
  geom_point()

n_fire_sims %>%
  ggplot(aes(z)) +
  geom_histogram()

# For each fire event, simulate fire size
n_fires <- sum(n_fire_sims$n_fire)

y_df <- Y %>%
  tbl_df %>%
  gather(sim_year, z) %>%
  mutate(pixel_idx = rep(rep(unique(st_covs$pixel_idx), times = 2), n_year),
         year = rep(unique(st_covs$year), each = nrow(st_covs) / n_year),
         dim = rep(rep(1:2, each = length(unique(st_covs$pixel_idx))), n_year),
         dim = factor(dim)) %>%
  full_join(st_covs)



fire_size_df <- tibble(pixel_idx = rep(n_fire_sims$pixel_idx, n_fire_sims$n_fire),
                       year = rep(n_fire_sims$year, n_fire_sims$n_fire))
fire_size_df <- y_df %>%
  filter(dim == 2) %>%
  dplyr::select(pixel_idx, year, z, dim) %>%
  right_join(fire_size_df)

stopifnot(nrow(fire_size_df) == n_fires)

# simulate log fire sizes
fire_size_df <- fire_size_df %>%
  mutate(z_obs = rnorm(n(), mean = z, sd = .1))

fire_size_df %>%
  ggplot(aes(z, z_obs)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_equal()

# Now find indices for each count and each fire size to match to mean matrix
count_idx <- rep(NA, nrow(n_fire_sims))
for (i in seq_along(count_idx)) {
  count_idx[i] <- which(master_idx_df$pixel_idx == n_fire_sims$pixel_idx[i] &
          master_idx_df$year == n_fire_sims$year[i] &
          master_idx_df$dim == 1)
}

size_idx <- rep(NA, nrow(fire_size_df))
for (i in seq_along(size_idx)) {
  size_idx[i] <-  which(master_idx_df$pixel_idx == fire_size_df$pixel_idx[i] &
                          master_idx_df$year == fire_size_df$year[i] &
                          master_idx_df$dim == 2)
}


stan_d <- list(n = length(unique(master_idx_df$pixel_idx)),
               T = n_year,
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

m_fit <- sampling(m_init, data = stan_d, cores = 3, chains = 3, iter = 400)
traceplot(m_fit, inc_warmup = TRUE)

traceplot(m_fit, pars = "sigma_size")
traceplot(m_fit, pars = "sigma_eta")

post <- rstan::extract(m_fit)

# Evaluate recovery of spatiotemporal basis coefficients
eta_df <- eta %>%
  reshape2::melt(varnames = c("year", "basis_dim"),
                 value.name = "true_eta")


post$eta %>%
  reshape2::melt(varnames = c("iter", "basis_dim", "year")) %>%
  tbl_df %>%
  group_by(basis_dim, year) %>%
  summarize(med = median(value),
            lo = quantile(value, .025),
            hi = quantile(value, .975)) %>%
  ungroup() %>%
  full_join(eta_df) %>%
  ggplot(aes(true_eta, med)) +
  geom_point() +
  facet_wrap(~basis_dim) +
  geom_segment(aes(xend = true_eta, y = lo, yend = hi)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  xlab("True basis coefficient") +
  ylab("Estimated basis coefficient") +
  coord_equal() +
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
  facet_wrap(~ idx, scales = "free") +
  geom_vline(aes(xintercept = true_value), linetype = "dashed")


post$mu %>%
  reshape2::melt(varnames = c("iter", "idx")) %>%
  mutate(pixel_idx = master_idx_df$pixel_idx[idx],
         year = master_idx_df$year[idx],
         dim = master_idx_df$dim[idx]) %>%
  full_join(y_df) %>%
  tbl_df %>%
  group_by(pixel_idx, year, dim) %>%
  summarize(med = median(value),
            lo = quantile(value, .025),
            hi = quantile(value, .975),
            true = unique(z)) %>%
  ggplot(aes(true, med)) +
  geom_segment(aes(xend = true, y = lo, yend = hi)) +
  facet_grid(dim ~ year) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")





# Visualize number of neighbors --------------------------------------
# theme_map <- theme(axis.line = element_blank(),
#                    axis.text.x = element_blank(),
#                    axis.text.y = element_blank(),
#                    axis.ticks = element_blank(),
#                    axis.title.x = element_blank(),
#                    axis.title.y = element_blank(),
#                    panel.background = element_blank(),
#                    panel.border = element_blank(),
#                    panel.grid.major = element_blank(),
#                    panel.grid.minor = element_blank(),
#                    plot.background = element_blank())

# ggplot(er_df, aes(long, lat, group = group)) +
#   geom_polygon(aes(fill = n_neighbors), color = NA) +
#   coord_equal() +
#   scale_fill_viridis() +
#   theme_map
#
#
# # Visualize the fire density in each ecoregion ----------------------------
# ggplot(er_df, aes(long, lat, group = group)) +
#   geom_polygon(aes(fill = log_fire_density), color = NA) +
#   coord_equal() +
#   labs(x = "Longitude", y = "Latitude") +
#   scale_fill_viridis("log(Fire density)") +
#   theme_map

# get count data (number of fires in each ecoregion X year)
data_summary <- d %>%
  group_by(us_l3name, fire_year, fire_mon) %>%
  summarize(n_fires = n()) %>%
  left_join(a_df) %>%
  full_join(all_ers) %>%
  ungroup() %>%
  dplyr::select(-area, -n_neighbors) %>%
  complete(us_l3name, fire_year, fire_mon,
           fill = list(n_fires = 0)) %>%
  full_join(all_ers) %>%
  mutate(cyear = c(scale(fire_year)),
         year = fire_year + 1 - min(fire_year),
         freg = factor(us_l3name,
                       levels = levels(factor(all_ers$us_l3name))),
         reg = as.numeric(freg),
         num_ym = as.numeric(factor(paste(fire_year,
                                          sprintf("%02d", fire_mon)))))

ymdf <- data_summary %>%
  distinct(fire_year, fire_mon, num_ym) %>%
  mutate(ymd = ymd(paste(fire_year,
                         sprintf("%02d", fire_mon),
                         sprintf("%02d", 1),
                         sep = "-")))


# get fire size data
fire_sizes <- d %>%
  dplyr::select(us_l3name, fire_year, fire_mon, fire_size) %>%
  mutate(cyear = c(scale(fire_year)),
         freg = factor(us_l3name,
                       levels = levels(data_summary$freg)),
         reg = as.numeric(freg),
         year = fire_year + 1 - min(fire_year)) %>%
  arrange(us_l3name, fire_year, fire_mon) %>%
  left_join(ymdf)


data_summary %>%
  ggplot(aes(as.numeric(fire_mon), n_fires / area, color = year)) +
  geom_path(alpha = .5) +
  scale_y_log10() +
  facet_wrap(~ us_l3name) +
  scale_color_viridis() +
  theme_minimal()
ggsave("fig/explore/fire-density-seasonality.pdf")

fire_sizes %>%
  ggplot(aes(as.numeric(fire_mon), fire_size, color = year)) +
  geom_jitter(alpha = .5, width = .3, size = .6) +
  scale_y_log10() +
  facet_wrap(~ us_l3name) +
  scale_color_viridis() +
  theme_minimal()
ggsave("fig/explore/fire-size-seasonality.pdf")

fire_sizes %>%
  ggplot(aes(ymd, fire_size, group = year)) +
  geom_point(alpha = .5) +
  scale_y_log10() +
  facet_wrap(~ us_l3name) +
  xlab("Year") +
  ylab("Fire size")
ggsave("fig/explore/fire-size-ts.pdf")


data_summary %>%
  left_join(ymdf) %>%
  ggplot(aes(ymd, n_fires/area)) +
  geom_line(alpha = .5) +
  scale_y_log10() +
  facet_wrap(~ us_l3name) +
  xlab("Year") +
  ylab("Fire density")
ggsave("fig/explore/fire-density-ts.pdf")
