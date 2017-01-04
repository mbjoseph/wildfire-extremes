source('R/clean_data.R')
library(ggplot2)
library(scales)
library(gridExtra)
library(spdep)
library(raster)
library(dplyr)
library(viridis)
library(tidyr)

d %>%
  ggplot(aes(x = discovery_, y = fire_size)) +
  geom_point() +
  scale_y_log10()

d %>%
  group_by(fire_year) %>%
  summarize(mean_size = mean(fire_size),
            median_size = median(fire_size),
            max_size = max(fire_size)) %>%
  gather(Measure, value, -fire_year) %>%
  ggplot(aes(x = fire_year, y = value, color = Measure)) +
  geom_point(aes(x = fire_year, y = fire_size), data = d, inherit.aes = FALSE) +
  geom_line() +
  theme_minimal() +
  scale_y_log10()


# Create spatial neighbors ------------------------------------------------
W <- poly2nb(ecoregions) %>%
  aggregate(ecoregions@data$US_L3NAME) %>%
  nb2mat(zero.policy = TRUE, style = 'B')

er_df <- d %>%
  group_by(us_l3name) %>%
  summarize(n_fires = n()) %>%
  right_join(er_df)


# get area & num_neighbors for each ecoregion and add to data frame
a_df <- data.frame(id = sapply(slot(ecoregions, "polygons"), slot, "ID"),
                   area = sapply(slot(ecoregions, "polygons"), slot, "area"),
                   US_L4NAME = ecoregions@data$US_L4NAME,
                   stringsAsFactors = FALSE) %>%
  group_by(US_L4NAME) %>%
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
  scale_fill_viridis() +
  theme_map

ggplot(er_df, aes(long, lat, group = group)) +
  geom_polygon(aes(fill = lfd), color = NA) +
  coord_equal() +
  labs(x = "Longitude", y = "Latitude") +
  scale_fill_viridis("log(Fire density)") +
  theme_map

all_ers <- er_df %>%
  group_by(us_l3name) %>%
  summarize(area = unique(area))
