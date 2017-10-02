
# Posterior predictive checks ---------------------------------------------

library(raster)
library(rstan)
library(sf)
library(rmapshaper)
library(tidyverse)
library(rstan)
library(modeest)
library(HDInterval)

source('R/02-explore.R')

ppred_df <- read_rds('ppred_df.rds') %>%
  mutate(total_burn_area = ifelse(is.na(total_burn_area),
                                  0,
                                  total_burn_area))

# compare empirical total burn areas to predicted
empirical_totals <- mtbs %>%
  as_tibble %>%
  group_by(NA_L3NAME, ym) %>%
  summarize(actual_total = sum(R_ACRES)) %>%
  ungroup %>%
  full_join(distinct(st_covs, NA_L3NAME, ym)) %>%
  complete(nesting(ym, NA_L3NAME),
           fill = list(actual_total = 0))

empirical_totals %>%
  ggplot(aes(ym, actual_total)) +
  geom_line() +
  facet_wrap(~ NA_L3NAME) +
  scale_y_log10()

predicted_totals <- ppred_df %>%
  group_by(idx) %>%
  summarize(pred_total = mlv(total_burn_area, method = 'venter')$M,
            median = median(total_burn_area),
            lo = hdi(total_burn_area, credMass = .9)[1],
            hi = hdi(total_burn_area, credMass = .9)[2])

burn_covs$idx <- 1:nrow(burn_covs)

predicted_totals <- predicted_totals %>%
  full_join(burn_covs)

total_ppred_df <- full_join(empirical_totals, predicted_totals)

total_ppred_df %>%
  mutate(Data = ifelse(ym < paste('Jan', cutoff_year),
                       'Train', 'Test')) %>%
  ggplot(aes(actual_total, pred_total, color = Data)) +
  theme_minimal() +
  geom_point(shape = 19, alpha = .5) +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(slope = 1, intercept = 0,
              linetype = 'dashed') +
  geom_segment(aes(xend = actual_total,
                   y = lo,
                   yend = hi),
               alpha = .05) +
  scale_color_gdocs() +
  facet_grid(Data ~ NA_L1NAME) +
  xlab('Empirical total burn area') +
  ylab('Predicted total burn area') +
  theme(legend.position = 'none') +
  coord_equal()
ggsave(filename = 'fig/ppc-burn-area.pdf', width = 20, height = 6)


# compare empirical maxima to predicted
empirical_maxima <- mtbs %>%
  as_tibble %>%
  group_by(NA_L3NAME, ym) %>%
  summarize(actual_max = max(R_ACRES)) %>%
  ungroup %>%
  full_join(distinct(st_covs, NA_L3NAME, ym)) %>%
  complete(nesting(ym, NA_L3NAME),
           fill = list(actual_max = 0))

empirical_maxima %>%
  ggplot(aes(ym, actual_max)) +
  geom_line() +
  facet_wrap(~ NA_L3NAME) +
  scale_y_log10()

predicted_maxima <- ppred_df %>%
  group_by(idx) %>%
  filter(!is.na(max_burn_area)) %>%
  summarize(pred_max = mlv(max_burn_area, method = 'venter')$M,
            median = median(max_burn_area),
            lo = hdi(max_burn_area, credMass = .9)[1],
            hi = hdi(max_burn_area, credMass = .9)[2]) %>%
  full_join(burn_covs)

max_ppred_df <- full_join(empirical_maxima, predicted_maxima) %>%
  filter(actual_max > 0)

max_ppred_df %>%
  mutate(Data = ifelse(ym < paste('Jan', cutoff_year),
                       'Train', 'Test')) %>%
  ggplot(aes(actual_max, pred_max, color = Data)) +
  theme_minimal() +
  scale_x_log10() +
  scale_y_log10() +
  geom_point() +
  geom_segment(aes(xend = actual_max,
                   y = lo,
                   yend = hi),
               alpha = .1) +
  geom_abline(slope = 1, intercept = 0,
              linetype = 'dashed') +
  scale_color_gdocs() +
  facet_grid(Data ~ NA_L1NAME) +
  xlab('Empirical max burn area') +
  ylab('Predicted max burn area') +
  theme(legend.position = 'none') +
  coord_equal()
ggsave(filename = 'fig/ppc-burn-max.pdf', width = 20, height = 6)
