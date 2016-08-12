source('R/thresholds/clean_data.R')
library(ggplot2)
library(ggmap)


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
grid.arrange(p1, p2)

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
  select(fire_size) %>%
  filter(fire_size > quantiles[which_quantile]) %>%
  unlist() %>%
  hist(main = paste('Distribution of values above the',
                     p[which_quantile],
                     'quantile'),
       breaks = 50)

extremes <- d %>%
  filter(fire_size > quantile(fire_size, .95))

hist(log(d$fire_size), breaks = 100)

bm <- get_map(location = 'Idaho, USA', zoom = 6, maptype = 'toner-background')

ggmap(bm) +
  geom_point(data = d, aes(x = longitude, y = latitude),
             alpha = .1, col = 'red') +
  stat_density2d(aes(x = longitude, y = latitude), data = d, h = 1)

ggmap(bm) +
  geom_point(data = extremes, aes(x = longitude, y = latitude),
             alpha = .5, col = 'red')

extremes %>%
  group_by(us_l4code) %>%
  summarize(n = n())
