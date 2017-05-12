library(scales)
library(gridExtra)
library(spdep)
library(viridis)
library(tidyverse)
library(lubridate)

source("R/01-clean_data.R")

# Visualize yearly fire size distributions ----------------------------
ymd <- d %>%
  group_by(fire_year, fire_mon) %>%
  summarize(mean_size = mean(fire_size),
            median_size = median(fire_size),
            max_size = max(fire_size)) %>%
  ungroup %>%
  tidyr::complete(fire_year, fire_mon,
           fill = list(mean_size = NA,
                       median_size = NA,
                       max_size = NA)) %>%
  mutate(yearmonth = paste(fire_year,
                           sprintf("%02d", fire_mon),
                           sep = "_")) %>%
  dplyr::select(-fire_year, -fire_mon) %>%
  mutate(num_ym = as.numeric(factor(yearmonth))) %>%
  gather(Measure, value, -yearmonth, -num_ym) %>%
  separate(yearmonth, c("year", "month"))


ymd %>%
  ggplot(aes(x = num_ym, y = value)) +
  geom_path() +
  facet_wrap(~ Measure) +
  scale_y_log10()

ymd %>%
  ggplot(aes(x = month, y = value, group = year)) +
  geom_path() +
  facet_wrap(~ Measure) +
  scale_y_log10()


# Create spatial neighbors ------------------------------------------------
W <- poly2nb(ecoregions) %>%
  aggregate(ecoregions@data$US_L3NAME) %>%
  nb2mat(zero.policy = TRUE, style = 'B')

er_df <- d %>%
  group_by(us_l3name) %>%
  summarize(n_fires = n()) %>%
  right_join(er_df)


## Add area & num_neighbors to data frame -------------------------------
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
         log_fire_density = log(n_fires) - log(area))

# get area of each ecoregion
all_ers <- er_df %>%
  group_by(us_l3name) %>%
  summarize(area = unique(area))



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
  distinct(fire_year, fire_mon, num_ym)


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
  ggplot(aes(fire_mon, n_fires / area, group = year)) +
  geom_path(alpha = .5) +
  scale_y_log10() +
  facet_wrap(~ us_l3name)

fire_sizes %>%
  ggplot(aes(fire_mon, fire_size, group = year)) +
  geom_point(alpha = .5) +
  scale_y_log10() +
  facet_wrap(~ us_l3name)

fire_sizes %>%
  ggplot(aes(num_ym, fire_size, group = year)) +
  geom_point(alpha = .5) +
  scale_y_log10() +
  facet_wrap(~ us_l3name)
