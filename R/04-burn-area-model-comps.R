
# Create tables for LOO-IC for each fit -----------------------------------
library(loo)
library(tidyverse)
library(cowplot)
library(ggthemes)

source('R/02-explore.R')
source('R/make-stan-d.R')

# possibly fetch from Amazon
#system("aws s3 cp s3://earthlab-mjoseph .. --recursive --exclude '*' --include '*.rds'")

model_fits <- list.files(pattern = '*.fit.*\\.rds')

burn_area_fits <- grep(model_fits, pattern = 'ba_', value = TRUE)

# for each model fit, produce a vector of the holdout log likelihoods
holdout_ba_loglik <- list()
train_ba_loglik <- list()
train_ba_rep <- list()
holdout_ba_rep <- list()

for (i in seq_along(burn_area_fits)) {
  post <- rstan::extract(read_rds(burn_area_fits[i]))
  holdout_ba_loglik[[i]] <- post$holdout_loglik_b %>%
    reshape2::melt(varnames = c('iter', 'idx')) %>%
    as_tibble %>%
    group_by(iter) %>%
    summarize(value = sum(value)) %>%
    mutate(model = burn_area_fits[i])

  train_ba_loglik[[i]] <- post$loglik_f %>%
    reshape2::melt(varnames = c('iter', 'idx')) %>%
    as_tibble %>%
    group_by(iter) %>%
    summarize(value = sum(value)) %>%
    mutate(model = burn_area_fits[i])

  holdout_ba_rep[[i]] <- post$holdout_rep %>%
    reshape2::melt(varnames = c('iter', 'idx')) %>%
    as_tibble %>%
    mutate(model = burn_area_fits[i])

  train_ba_rep[[i]] <- post$size_rep %>%
    reshape2::melt(varnames = c('iter', 'idx')) %>%
    as_tibble %>%
    mutate(model = burn_area_fits[i])
}

holdout_ba_loglik <- bind_rows(holdout_ba_loglik) %>%
  mutate(train = FALSE)

train_ba_loglik <- bind_rows(train_ba_loglik) %>%
  mutate(train = TRUE)

# data frame for plotting colors
cols <- c('Gamma' = 'purple',
          'Weibull' = 'orange',
          'Generalized Pareto' = 'red',
          'Tapered Pareto' = 'limegreen',
          'Lognormal' = 'dodgerblue')


holdout_ba_loglik %>%
  full_join(train_ba_loglik) %>%
  mutate(train = ifelse(train == TRUE, 'train', 'test'),
         Distribution = case_when(
           grepl('gamma', .$model) ~ 'Gamma',
           grepl('lognormal', .$model) ~ 'Lognormal',
           grepl('_pareto', .$model) ~ 'Generalized Pareto',
           grepl('_tpareto', .$model) ~ 'Tapered Pareto',
           grepl('weibull', .$model) ~ 'Weibull'
         )) %>%
  spread(train, value) %>%
  ggplot(aes(x = train, y = test, color = Distribution)) +
  geom_point(alpha = .05) +
  xlab('Log likelihood: training set') +
  ylab('Log likelihood: test set') +
  scale_color_manual('Burn area distribution', values = cols) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0)) +
  ylim(c(-26000, -22200)) +
  theme_minimal()
ggsave(filename = 'fig/loglik.png', width = 6, height = 4)

## Generate plots to evaluate distributional assumptions
train_ba_rep <- train_ba_rep %>%
  bind_rows

# Whole distribution
den_plot <- train_ba_rep %>%
  filter(iter < 101) %>%
  mutate(Distribution = case_when(
    grepl('gamma', .$model) ~ 'Gamma',
    grepl('lognormal', .$model) ~ 'Lognormal',
    grepl('_pareto', .$model) ~ 'Generalized Pareto',
    grepl('_tpareto', .$model) ~ 'Tapered Pareto',
    grepl('weibull', .$model) ~ 'Weibull'
  )) %>%
  ggplot(aes(x = value,
             color = Distribution,
             group = interaction(Distribution, iter))) +
  theme_minimal() +
  stat_density(aes(color = Distribution), geom="line", alpha = .1, position = 'identity') +
  scale_x_log10() +
  scale_color_manual('Burn area distribution', values = cols) +
  facet_wrap(~Distribution, nrow = 1) +
  stat_density(geom = 'line', position = 'identity',
               inherit.aes = FALSE,
               aes(x = actual_sizes),
               data = tibble(actual_sizes = stan_d$sizes)) +
  theme(legend.position = 'none') +
  coord_cartesian(xlim = c(min(stan_d$sizes), 700000)) +
  geom_rug(inherit.aes = FALSE,
           aes(x = actual_sizes),
           data = tibble(actual_sizes = stan_d$sizes),
           alpha = .1) +
  xlab('Burn area exceedance (acres)') +
  ylab('Density')

# "Survival" plots for top 3 models
obs_ccdf <- tibble(values = sort(stan_d$sizes)) %>%
  mutate(ccdf = seq(n(), 1, -1) / n()) %>%
  filter(ccdf <= .01)
gc()

tail_plot <- train_ba_rep %>%
  filter(iter < 1001) %>%
  group_by(iter, model) %>%
  arrange(iter, model, value) %>%
  mutate(ccdf = seq(n(), 1, -1) / n()) %>%
  filter(ccdf <= .01) %>%
  group_by(ccdf, model) %>%
  summarize(lo = quantile(value, .025),
            q25 = quantile(value, .25),
            med = median(value),
            q75 = quantile(value, .75),
            hi = quantile(value, .975)) %>%
  ungroup %>%
  mutate(Distribution = case_when(
    grepl('gamma', .$model) ~ 'Gamma',
    grepl('lognormal', .$model) ~ 'Lognormal',
    grepl('_pareto', .$model) ~ 'Generalized Pareto',
    grepl('_tpareto', .$model) ~ 'Tapered Pareto',
    grepl('weibull', .$model) ~ 'Weibull'
  )) %>%
  ggplot(aes(med, ccdf, color = Distribution)) +
  theme_minimal() +
  scale_color_manual('Burn area distribution', values = cols) +
  geom_point(shape = 1) +
  geom_errorbarh(aes(xmin = lo, xmax = hi)) +
  geom_errorbarh(aes(xmin = q25, xmax = q75), size = 1.25) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~Distribution, nrow = 1) +
  coord_cartesian(xlim = c(1e5, 1e8),
                  ylim = c(1e-4, .01)) +
  theme(legend.position = 'none') +
  geom_point(data = obs_ccdf, color = 'black', aes(x = values), shape = 1, size = .5) +
  xlab('Burn area exceedance (acres)') +
  ylab('Complementary CDF') +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )


# Maxima and cumulative burn area -----------------------------------------
train_max <- train_ba_rep %>%
  group_by(iter, model) %>%
  summarize(max_pred = max(value),
            cum_pred = sum(value)) %>%
  mutate(train = TRUE)

holdout_max <- holdout_ba_rep %>%
  bind_rows  %>%
  group_by(iter, model) %>%
  summarize(max_pred = max(value),
            cum_pred = sum(value)) %>%
  mutate(train = FALSE)

gc()

actual_maxima <- tibble(train = max(stan_d$sizes),
                        test = max(stan_d$holdout_b))

actual_sums <- tibble(train = sum(stan_d$sizes),
                      test = sum(stan_d$holdout_b))

max_plot <- train_max %>%
  full_join(holdout_max) %>%
  dplyr::select(-cum_pred) %>%
  mutate(train = ifelse(train, 'train', 'test')) %>%
  spread(train, max_pred) %>%
  ungroup() %>%
  mutate(Distribution = case_when(
    grepl('gamma', .$model) ~ 'Gamma',
    grepl('lognormal', .$model) ~ 'Lognormal',
    grepl('_pareto', .$model) ~ 'Generalized Pareto',
    grepl('_tpareto', .$model) ~ 'Tapered Pareto',
    grepl('weibull', .$model) ~ 'Weibull'
  )) %>%
  ggplot(aes(train, test, color = Distribution)) +
  theme_minimal() +
  scale_color_manual('Burn area distribution', values = cols) +
  geom_point(alpha = .2, shape = 1) +
  xlab('max(acres): training data') +
  ylab('max(acres): test data') +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~ Distribution, nrow = 1) +
  theme(legend.position = 'none') +
  geom_point(data = actual_maxima, color = 'black') +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

sum_plot <- train_max %>%
  full_join(holdout_max) %>%
  dplyr::select(-max_pred) %>%
  mutate(train = ifelse(train, 'train', 'test')) %>%
  spread(train, cum_pred) %>%
  ungroup() %>%
  mutate(Distribution = case_when(
    grepl('gamma', .$model) ~ 'Gamma',
    grepl('lognormal', .$model) ~ 'Lognormal',
    grepl('_pareto', .$model) ~ 'Generalized Pareto',
    grepl('_tpareto', .$model) ~ 'Tapered Pareto',
    grepl('weibull', .$model) ~ 'Weibull'
  )) %>%
  ggplot(aes(train, test, color = Distribution)) +
  theme_minimal() +
  scale_color_manual('Burn area distribution', values = cols) +
  geom_point(alpha = .2, shape = 1) +
  xlab('sum(acres): training data') +
  ylab('sum(acres): test data') +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~ Distribution, nrow = 1) +
  theme(legend.position = 'none') +
  geom_point(data = actual_sums, color = 'black') +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  # add dashed lines to indicate area of contiguous U.S.
  geom_vline(xintercept = 1996726201.6, linetype = 'dotted') +
  geom_hline(yintercept = 1996726201.6, linetype = 'dotted')

plot_grid(den_plot, tail_plot, max_plot, sum_plot, nrow = 4, labels = 'auto')
ggsave(filename = 'fig/ppc-density-funs.png', width = 12, height = 10)
