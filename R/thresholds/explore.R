source('R/thresholds/clean_data.R')
library(ggplot2)
library(scales)
library(gridExtra)
library(spdep)
library(raster)
library(dplyr)

d %>%
  ggplot(aes(x = discovery_, y = fire_size)) +
  geom_point() +
  scale_y_log10()

summ_d <- d %>%
  group_by(fire_year) %>%
  summarize(mean_size = mean(fire_size),
            median_size = median(fire_size),
            n = n(),
            max_size = max(fire_size))

p1 <- summ_d %>%
  ggplot(aes(fire_year)) +
  geom_point(aes(x = fire_year, y = fire_size), data = d) +
  geom_line(aes(y = mean_size), col = 'blue') +
  geom_line(aes(y = median_size), col = 'red') +
  theme_minimal() +
  scale_y_log10()

p2 <- summ_d %>%
  arrange(fire_year) %>%
  ggplot(aes(n, max_size, label = fire_year)) +
  geom_path(arrow = arrow(), alpha = .2) +
  geom_text() +
  theme_minimal()


grid.arrange(p1, p2, nrow = 1)

# trying to find a threshold
plot(sort(d$fire_size))
p <- c(.9, .95, .99, .999)
quantiles <- quantile(d$fire_size, probs = p)
abline(h = quantiles, col = 2, lty = 2)
text(x = 1000, y = quantiles,
     labels = paste(100 * p, '% quantile'))

# make histogram for values above a threshold
which_quantile <- 1
d %>%
  dplyr::select(fire_size) %>%
  filter(fire_size > quantiles[which_quantile]) %>%
  unlist() %>%
  hist(main = paste('Distribution of values above the',
                     p[which_quantile],
                     'quantile'),
       breaks = 50)

extremes <- d %>%
  filter(fire_size > quantile(fire_size, .95))

hist(log(d$fire_size), breaks = 100)




# Create spatial neighbors ------------------------------------------------
coords <- coordinates(ecoregions)
Wpoly <- poly2nb(ecoregions)
Wagg <- aggregate(Wpoly, ecoregions@data$US_L3NAME)
W <- nb2mat(Wagg, zero.policy = TRUE, style = 'B')

n_fires <- d %>%
  group_by(us_l3name) %>%
  summarize(n_fires = n())

er_df <- left_join(er_df, n_fires)

# get area & num_neighbors for each ecoregion and add to data frame
a_df <- data.frame(id = sapply(slot(ecoregions, "polygons"), slot, "ID"),
                   area = sapply(slot(ecoregions, "polygons"), slot, "area"),
                   US_L3NAME = ecoregions@data$US_L3NAME,
                   stringsAsFactors = FALSE) %>%
  group_by(US_L3NAME) %>%
  summarize(area = sum(area))
a_df$n_neighbors <- rowSums(W)
names(a_df) <- tolower(names(a_df))

er_df <- left_join(er_df, a_df) %>%
  mutate(fire_den = n_fires / area,
         lfd = log(n_fires) - log(area))

theme_map <- theme(axis.line = element_blank(),
                   axis.text.x = element_blank(),
                   axis.text.y = element_blank(),
                   axis.ticks = element_blank(),
                   axis.title.x = element_blank(),
                   axis.title.y = element_blank(),
                   panel.background = element_blank(),
                   panel.border = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   plot.background = element_blank())

ggplot(er_df, aes(long, lat, group = group)) +
  geom_polygon(aes(fill = n_neighbors)) +
  geom_polygon(color = 'black', alpha = .1, size = .1) +
  coord_equal() +
  scale_fill_gradientn(colors = c('black', 'darkblue', 'blue', 'dodgerblue',
                                  'lightblue')) +
  theme_map

ggplot(er_df, aes(long, lat, group = group)) +
  geom_polygon(aes(fill = lfd), color = NA) +
  coord_equal() +
  labs(x = "Longitude", y = "Latitude") +
  scale_fill_gradientn(colors = c('black', muted('red'), 'red'),
                       na.value = 'black') +
  theme_map

er_df %>%
  group_by(group) %>%
  summarize(fire_den = unique(fire_den)) %>%
  ggplot(aes(x = fire_den)) +
  geom_histogram(bins = 20)

# simulate realizations from CAR prior
D <- diag(rowSums(W))
invD <- solve(D)
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
isqrtD <- solve(expm::sqrtm(D))
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
