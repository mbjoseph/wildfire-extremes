library(tidyverse)
library(patchwork)
library(grid)
st_covs <- read_rds('data/processed/st_covs.rds')


# Load and parse burn area partial effs -----------------------------------
ba_partials <- 'data/processed/burn_area_partial_all.rds' %>%
  read_rds

for (i in seq_along(ba_partials)) {
  out <- ba_partials[[i]] %>%
    bind_rows %>%
    left_join(st_covs) %>%
    mutate(NA_L1NAME = tolower(NA_L1NAME),
           NA_L1NAME = factor(tools::toTitleCase(NA_L1NAME)),
           NA_L1NAME = fct_reorder(NA_L1NAME, rmin))
  ba_partials[[i]] <- out
}
ba_partials <- bind_rows(ba_partials, .id = 'covariate')

# Load and parse count partial effs ---------------------------------------
c_partials <- read_csv('data/processed/count_partial_effect_plot_data.csv') %>%
  mutate(response = 'Fire counts') %>%
  select(covariate_value, med, NA_L3NAME, fancy_name, response, covariate)


covs_to_plot <- c('tmmx', 'rmin')


ba_df <- ba_partials %>%
  filter(covariate %in% covs_to_plot) %>%
  mutate(covariate_value = ifelse(covariate == 'tmmx', tmmx, rmin), 
         response = 'Burned area') %>%
  select(covariate_value, med, NA_L3NAME, response, covariate) %>%
  left_join(distinct(c_partials, covariate, fancy_name))



# Visualize ---------------------------------------------------------------
title_size <- 11
margin_size <- unit(c(rep(.01, 4)), "cm")
p1 <- ba_df %>%
  ggplot(aes(covariate_value, y = med)) +
  theme_minimal() +
  geom_line(aes(group = NA_L3NAME), alpha = .2) +
  ylab('Partial effect') +
  facet_wrap(~fancy_name, scales = 'free_x', strip.position = 'bottom', nrow = 1) +
  xlab('') +
  theme(panel.grid.minor = element_blank(),
        legend.position = 'bottom') + 
  scale_y_continuous(limits = c(-.1, 1.5)) + 
  ggtitle('Wildfire Burn Area') + 
  theme(plot.title = element_text(size=title_size), 
        plot.margin = margin_size)

p2 <- c_partials %>%
  filter(covariate %in% covs_to_plot) %>%
  ggplot(aes(covariate_value, y = med)) +
  theme_minimal() +
  geom_line(aes(group = NA_L3NAME), alpha = .2) +
  ylab('Partial effect') +
  facet_wrap(~fancy_name, scales = 'free_x', strip.position = 'bottom', nrow = 1) +
  xlab('') +
  theme(panel.grid.minor = element_blank(),
        legend.position = 'bottom') + 
  scale_y_continuous(breaks = c(-3, 0, 3)) + 
  ggtitle("Wildfire Count") + 
  theme(plot.title = element_text(size=title_size), 
        plot.margin = margin_size)

p1 / p2
