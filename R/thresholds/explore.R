source('R/thresholds/clean_data.R')
library(ggplot2)
library(scales)
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

library(gridExtra)
grid.arrange(p1, p2, nrow = 1)

d %>%
  ggplot(aes(x = fire_year, y = fire_size)) +
  geom_point() +
  facet_wrap(~stat_cau_1) +
  stat_smooth() +
  scale_y_log10()

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

bm <- get_map(location = 'Kansas, USA', zoom = 4, maptype = 'toner-background')

ggmap(bm) +
  geom_point(data = d, aes(x = longitude, y = latitude),
             alpha = .1, col = 'red') +
  stat_density2d(aes(x = longitude, y = latitude), data = d, h = 1)

ggmap(bm) +
  geom_point(data = extremes, aes(x = longitude, y = latitude),
             alpha = .5, col = 'red')

library(spdep)
library(raster)
coords <- coordinates(ecoregions)
nb <- knn2nb(knearneigh(coords, k = 1))
max_d <- max(unlist(nbdists(nb, coords)))
nb <- dnearneigh(coords, 0, max_d)
W <- nb2mat(nb, zero.policy = TRUE)

plot(ecoregions, border = 'grey')
plot(nb, coords, add = TRUE, points = FALSE,
     col = alpha('blue', .1), lwd = .5)
points(spd, cex = .2, col = alpha('red', .5))

n_fires <- d %>%
  group_by(us_l4name) %>%
  summarize(n_fires = n()) %>%
  mutate(US_L4NAME = us_l4name)
er_df <- left_join(er_df, n_fires)

a_df <- data.frame(id = sapply(slot(ecoregions, "polygons"), slot, "ID"),
                   area = sapply(slot(ecoregions, "polygons"), slot, "area"),
                   stringsAsFactors = FALSE)
er_df <- left_join(er_df, a_df) %>%
  mutate(fire_den = n_fires / area,
         lfd = log(n_fires) - log(area))

ggplot(er_df, aes(long, lat, group = group)) +
  geom_polygon(aes(fill = US_L4NAME)) +
  geom_polygon(color = 'black', alpha = .1, size = .1) +
  coord_equal() +
  labs(x = "Longitude", y = "Latitude") +
  theme(legend.position = 'none')

ggplot(er_df, aes(long, lat, group = group)) +
  geom_polygon(aes(fill = lfd), color = NA) +
  coord_equal() +
  labs(x = "Longitude", y = "Latitude") +
  scale_fill_gradient2(low = 'black', mid = 'white', high = 'red',
                       na.value = 'black') +
  theme_minimal()

er_df %>%
  group_by(id) %>%
  summarize(fire_den = unique(fire_den)) %>%
  ggplot(aes(x = log(fire_den))) +
  geom_histogram(bins = 100)

# simulate realizations from CAR prior
rho <- .8
n <- nrow(W)
n_islands <- sum(rowSums(W) == 0)
tau <- 2
Sigma <- solve(diag(n) - rho * W) %*%
  diag(tau ^ 2 / rowSums(W))
plot(raster(Sigma[1:100, 1:100]))
L <- t(chol(Sigma))

sim_df <- data.frame(id = rownames(W),
                     ysim = L %*% rnorm(n))

er_df <- full_join(er_df, sim_df)

ggplot(er_df, aes(long, lat, group = group)) +
  geom_polygon(aes(fill = ysim), color = NA) +
  coord_equal() +
  labs(x = "Longitude", y = "Latitude") +
  scale_fill_gradient2(low = 'black', mid = 'red', high = 'yellow') +
  theme_minimal()
