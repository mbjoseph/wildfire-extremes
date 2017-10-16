
# Posterior predictive checks ---------------------------------------------

library(raster)
library(rstan)
library(sf)
library(rmapshaper)
library(tidyverse)
library(rstan)

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

st_covs$idx <- 1:nrow(st_covs)

# compute Bayesian posterior predictive p-values
pvdf <- empirical_totals %>%
  full_join(select(st_covs, idx, ym, NA_L3NAME)) %>%
  full_join(select(ppred_df, idx, total_burn_area)) %>%
  filter(actual_total > 0) %>%
  group_by(ym, NA_L3NAME) %>%
  summarize(pval = mean(total_burn_area >= actual_total)) %>%
  mutate(Data = ifelse(ym < paste('Jan', cutoff_year),
                       'Train', 'Test'),
         extreme = pval < 0.05 | pval > 0.95)
gc()
hist(pvdf$pval, breaks = 100)

pvdf %>%
  ggplot(aes(ym, pval, color = Data,
             alpha = extreme)) +
  geom_point() +
  facet_wrap(~ NA_L3NAME) +
  geom_hline(yintercept = .5, linetype = 'dashed') +
  geom_hline(yintercept = c(.05, .95), linetype = 'dotted', alpha = .4) +
  scale_color_gdocs()

# graphical predictive check
predicted_totals <- ppred_df %>%
  group_by(idx) %>%
  summarize(pred_total = median(total_burn_area),
            lo = quantile(total_burn_area, .05),
            hi = quantile(total_burn_area, .95))


predicted_totals <- predicted_totals %>%
  full_join(st_covs)

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


## compare total burn area over entire record --------------------
pred_gt <- ppred_df %>%
  group_by(iter) %>%
  summarize(grand_total_pred = sum(total_burn_area)) %>%
  ungroup

actual_gt <- sum(mtbs$R_ACRES)

pred_gt %>%
  ggplot(aes(log(grand_total_pred))) +
  geom_histogram(bins = 300) +
  geom_vline(xintercept = log(actual_gt), color = 'red', linetype = 'dashed') +
  scale_y_log10()

mean(pred_gt$grand_total_pred >= actual_gt)
# oh farts, we are predicting way more total burn area than what's observed...


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
  filter(!is.na(max_burn_area)) %>%
  group_by(idx) %>%
  summarize(pred_max = median(max_burn_area),
            lo = quantile(max_burn_area, .05),
            hi = quantile(max_burn_area, .95)) %>%
  full_join(st_covs)


max_ppred_df <- full_join(empirical_maxima, predicted_maxima) %>%
  filter(actual_max > 0)

max_ppred_df %>%
  mutate(Data = ifelse(ym < paste('Jan', cutoff_year),
                       'Train', 'Test')) %>%
  ggplot(aes(actual_max, pred_max, color = Data)) +
  theme_minimal() +
  scale_x_log10() +
  scale_y_log10() +
  geom_point(shape = 1) +
  geom_segment(aes(xend = actual_max,
                   y = lo,
                   yend = hi),
               alpha = .2) +
  geom_abline(slope = 1, intercept = 0,
              linetype = 'dashed') +
  scale_color_gdocs() +
  facet_grid(Data ~ NA_L1NAME) +
  xlab('Empirical max burn area') +
  ylab('Predicted max burn area') +
  theme(legend.position = 'none') +
  coord_equal()
ggsave(filename = 'fig/ppc-burn-max.pdf', width = 20, height = 6)



# Fraction of the year with fires -----------------------------------------
empirical_pwf <- count_df %>%
  group_by(FIRE_YEAR, NA_L3NAME) %>%
  summarize(pwf = mean(n_fire > 0)) %>%
  ungroup %>%
  rename(year = FIRE_YEAR)

empirical_pwf %>%
  ggplot(aes(year, pwf)) +
  facet_wrap(~NA_L3NAME) +
  geom_line()

predicted_pwf <- ppred_df %>%
  full_join(st_covs) %>%
  group_by(iter, year, NA_L3NAME) %>%
  summarize(pwf = mean(n_events > 0)) %>%
  ungroup %>%
  group_by(year, NA_L3NAME) %>%
  summarize(pred_pwf = median(pwf),
            lo = quantile(pwf, .025),
            hi = quantile(pwf, .975)) %>%
  ungroup

full_join(empirical_pwf, predicted_pwf) %>%
  mutate(is_train = year < cutoff_year) %>%
  filter(!is_train) %>%
  ggplot(aes(year)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = is_train),
              alpha = .5) +
  geom_point(aes(y = pwf)) +
  facet_wrap(~NA_L3NAME) +
  scale_fill_gdocs() +
  xlab('Year') +
  ylab('Proportion of year with fires')
ggsave(filename = 'fig/ppc-p-year.pdf', width = 20, height = 9)

# Investigate why the model is overpredicting the fraction of year
predicted_counts <- ppred_df %>%
  full_join(select(st_covs, idx, NA_L3NAME, ym)) %>%
  group_by(NA_L3NAME, ym) %>%
  summarize(pred_num = median(n_events),
            lo = quantile(n_events, .05),
            hi = quantile(n_events, .95)) %>%
  ungroup

predicted_counts %>%
  full_join(count_df) %>%
  mutate(is_train = year < cutoff_year) %>%
  filter(!is_train) %>%
  ggplot(aes(ym)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = is_train),
              alpha = .5) +
  geom_point(aes(y = n_fire), size = .5) +
  facet_wrap(~NA_L3NAME) +
  scale_fill_gdocs() +
  xlab('Year') +
  ylab('Number of fires') +
  scale_y_log10()

