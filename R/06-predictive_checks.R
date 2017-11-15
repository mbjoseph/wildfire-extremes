
# Posterior predictive checks ---------------------------------------------

library(raster)
library(rstan)
library(sf)
library(rmapshaper)
library(tidyverse)
library(rstan)

source('R/02-explore.R')

ba_df <- read_rds('ppc_ba.rds') %>%
  filter(!is.infinite(value)) %>%
  mutate(iteration = parse_number(iteration))


ba_df %>%
  ggplot(aes(x = log(value - 1e3))) +
  scale_x_log10(limits = c(2,  20)) +
  geom_line(aes(group = iteration), alpha = .2, color = alpha(1, .1),
               stat = 'density') +
  geom_line(aes(x = log(R_ACRES - 1e3)),
               data = mtbs, color = 'dodgerblue', size = 2, alpha = .8,
               stat = 'density') +
  geom_jitter(aes(x = log(R_ACRES - 1e3), y = -.1), width = 0, height = .1, data = mtbs,
              color = 'dodgerblue', alpha = .5, shape = 1) +
  xlab('Logarithm of exceedence') +
  ylab('Density')

log_maxima <- ba_df %>%
  group_by(iteration) %>%
  summarize(max_log_burn_area = max(log(value)))

log_maxima %>%
  ggplot(aes(max_log_burn_area)) +
  geom_line(stat = 'density') +
  scale_x_log10() +
  geom_vline(xintercept = max(log(mtbs$R_ACRES)), col = 'dodgerblue')

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

# graphical predictive check
predicted_totals <- ppred_df %>%
  group_by(idx) %>%
  summarize(pred_total = median(total_burn_area),
            lo = quantile(total_burn_area, .05),
            hi = quantile(total_burn_area, .95),
            q25 = quantile(total_burn_area, .25),
            q75 = quantile(total_burn_area, .75))


predicted_totals <- predicted_totals %>%
  full_join(st_covs)

total_ppred_df <- full_join(empirical_totals, predicted_totals)

total_ppred_df %>%
  mutate(Data = ifelse(ym < paste('Jan', cutoff_year),
                       'Train', 'Test')) %>%
  filter(actual_total > 0, pred_total > 0) %>%
  ggplot(aes(actual_total, pred_total, color = Data)) +
  theme_minimal() +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(slope = 1, intercept = 0,
              linetype = 'dashed') +
  geom_segment(aes(xend = actual_total,
                   y = lo,
                   yend = hi),
               alpha = .1) +
  geom_segment(aes(xend = actual_total,
                   y = q25,
                   yend = q75),
               alpha = .2, size = 1.1) +
  geom_point(shape = 19, alpha = .5) +
  scale_color_gdocs() +
  facet_wrap(~ Data) +
  xlab('Empirical total burn area') +
  ylab('Predicted total burn area') +
  theme(legend.position = 'none')
ggsave(filename = 'fig/ppc-burn-area.pdf', width = 26, height = 10)

# check interval coverage (~ 98% coverage for a 90% interval! :)
total_ppred_df %>%
  mutate(Data = ifelse(ym < paste('Jan', cutoff_year),
                       'Train', 'Test')) %>%
  mutate(lo_lt_real = lo <= actual_total,
         hi_gt_real = hi >= actual_total,
         actual_in_interval = lo_lt_real & hi_gt_real) %>%
  group_by(Data) %>%
  summarize(coverage = mean(actual_in_interval))

## compare total burn area over entire record --------------------
pred_gt <- ppred_df %>%
  group_by(iter) %>%
  summarize(grand_total_pred = sum(total_burn_area)) %>%
  ungroup

actual_gt <- sum(mtbs$R_ACRES)

pred_gt %>%
  ggplot(aes(log(grand_total_pred))) +
  geom_density() +
  geom_vline(xintercept = log(actual_gt), color = 'red', linetype = 'dashed') +
  scale_x_log10() +
  annotation_logticks(sides = 'b') +
  xlab('log(Grand total burn area)') +
  ylab('Posterior density')

mean(pred_gt$grand_total_pred >= actual_gt)

gc()


# What about the total number of fires ------------------------------------
pred_n <- ppred_df %>%
  group_by(iter) %>%
  summarize(grand_total_n = sum(n_events)) %>%
  ungroup

actual_n <- nrow(mtbs)

pred_n %>%
  ggplot(aes(log(grand_total_n))) +
  geom_density() +
  geom_vline(xintercept = log(actual_n), color = 'red', linetype = 'dashed') +
  xlab('log(Grand total burn number of fires)') +
  ylab('Posterior density')
mean(pred_gt$grand_total_pred >= actual_gt)



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
ggsave(filename = 'fig/ppc-burn-max.pdf', width = 26, height = 10)

# check interval coverage (90% train, 87% test for 90% interval :)
max_ppred_df %>%
  mutate(Data = ifelse(ym < paste('Jan', cutoff_year),
                       'Train', 'Test')) %>%
  na.omit %>%
  mutate(lo_lt_real = lo <= actual_max,
         hi_gt_real = hi >= actual_max,
         actual_in_interval = lo_lt_real & hi_gt_real) %>%
  group_by(Data) %>%
  summarize(coverage = mean(actual_in_interval))
gc()






# Fraction of the year with fires -----------------------------------------
rm(list = ls())
# restart
gc()

library(rstan)
library(sf)
library(rmapshaper)
library(tidyverse)
library(rstan)

source('R/02-explore.R')
st_covs$idx <- 1:nrow(st_covs)

ppred_df <- read_rds('ppred_df.rds') %>%
  mutate(total_burn_area = ifelse(is.na(total_burn_area),
                                  0,
                                  total_burn_area))

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
  full_join(dplyr::select(st_covs, idx, ym, year, NA_L3NAME)) %>%
  mutate(Data = ifelse(ym < paste('Jan', cutoff_year),
                       'Train', 'Test')) %>%
  filter(Data == 'Test') %>%
  group_by(iter, year, NA_L3NAME) %>%
  summarize(pwf = mean(n_events > 0)) %>%
  ungroup %>%
  group_by(year, NA_L3NAME) %>%
  summarize(pred_pwf = median(pwf),
            lo = quantile(pwf, .025),
            hi = quantile(pwf, .975)) %>%
  ungroup

right_join(empirical_pwf, predicted_pwf) %>%
  mutate(is_train = year < cutoff_year) %>%
  filter(!is_train) %>%
  ggplot(aes(year)) +
  geom_ribbon(aes(ymin = lo, ymax = hi),
              alpha = .5) +
  geom_point(aes(y = pwf)) +
  facet_wrap(~NA_L3NAME) +
  scale_fill_gdocs() +
  xlab('Year') +
  ylab('Proportion of year with fires') +
  ggtitle('Posterior predictive check for the proportion of year with 1 or more fires, test data')
ggsave(filename = 'fig/ppc-p-year.pdf', width = 20, height = 9)

# Investigate why the model is overpredicting the fraction of year
predicted_counts <- ppred_df %>%
  full_join(dplyr::select(st_covs, idx, NA_L3NAME, ym)) %>%
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
  geom_ribbon(aes(ymin = lo, ymax = hi),
              alpha = .5) +
  geom_point(aes(y = n_fire), size = .5) +
  facet_wrap(~NA_L3NAME) +
  scale_fill_gdocs() +
  xlab('Year') +
  ylab('Number of fires') +
  scale_y_log10() +
  theme_bw() +
  ggtitle('Posterior predictive check for the number of fires, test data')
ggsave(filename = 'fig/ppc-n-fire.pdf', width = 20, height = 9)

# proportion of empirical points in prediction interval
predicted_counts %>%
  full_join(count_df) %>%
  mutate(is_train = year < cutoff_year) %>%
  filter(!is_train) %>%
  group_by(NA_L3NAME) %>%
  summarize(p_in = mean(n_fire >= lo & n_fire <= hi)) %>%
  ggplot(aes(reorder(NA_L3NAME, p_in), p_in)) +
  geom_point() +
  coord_flip() +
  ylab('Interval converage: number of fires') +
  xlab('Ecoregion')



# Cumulative burn area check ----------------------------------------------
mtbs_hist <- read_rds('mtbs_hist.rds')
hist_df <- read_rds('ppc_hist.rds')


mtbs_hist_df <- tibble(count = mtbs_hist$counts,
                       midpoint = mtbs_hist$mids) %>%
  mutate(approx_burn_area = count * midpoint,
         cum_burn_area = cumsum(approx_burn_area),
         iter = 1,
         p_burn_area = cum_burn_area / max(cum_burn_area))


hist_df %>%
  ggplot(aes(log(midpoint), log(cum_burn_area), group = iter)) +
  geom_line(alpha = .1) +
  geom_line(data = mtbs_hist_df, color = 'red') +
  coord_cartesian(ylim = c(10.5, 21),
                  xlim = c(5, 16)) +
  xlab('log(Burn area)') +
  ylab('log(Cumulative burn area)')

# normalized by total burn area
hist_df %>%
  group_by(iter) %>%
  mutate(total_burn_area = max(cum_burn_area),
         p_burn_area = cum_burn_area / total_burn_area) %>%
  ggplot(aes(log(midpoint), p_burn_area, group = iter)) +
  geom_line(alpha = .1) +
  geom_line(data = mtbs_hist_df, color = 'red') +
  xlab('log(Burn area)') +
  scale_x_log10() +
  ylab('Cumulative proportion of total burn area') +
  theme_bw()




# Compare the empirical vs. predicted CDF for # fires -----------------------------
max_breaks <- 200
count_hist <- hist(count_df$n_fire, breaks = max_breaks, right = FALSE)
count_hist_df <- tibble(count = count_hist$counts,
                        midpoint = count_hist$mids,
                        p_zero = mean(count_df$n_fire == 0)) %>%
  mutate(cum_count = cumsum(count),
         iter = 1,
         p_count = cum_count / max(cum_count))

ggplot(count_hist_df, aes(midpoint, p_count)) +
  geom_line() +
  geom_point(x = 0, aes(y = p_zero))

# generate the histograms for each iteration of the posterior and compare
pred_hist_df <- ppred_df %>%
  split(.$iter) %>%
  map(~ hist(.x$n_events, breaks = max_breaks, right = FALSE, plot = FALSE)) %>%
  map(~ tibble(midpoint = .x$mids,
               count = .x$counts) %>%
        mutate(cum_count = cumsum(count),
               p_count = cum_count / max(cum_count))) %>%
  bind_rows(.id = 'iter')

p_zero_df <- ppred_df %>%
  group_by(iter) %>%
  summarize(p_zero = mean(n_events == 0))

pred_hist_df %>%
  ggplot(aes(midpoint, p_count, group = iter)) +
  geom_line(alpha = .05) +
  geom_line(data = count_hist_df, color = 'red') +
  geom_point(aes(x = 0, y = p_zero), alpha = .05,
             data = p_zero_df) +
  geom_point(aes(x = 0, y = p_zero), data = count_hist_df, col ='red') +
  coord_cartesian(xlim = c(0, 20)) +
  xlab('Number of fires') +
  ylab('CDF(Number of fires)')
