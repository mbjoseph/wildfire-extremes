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
