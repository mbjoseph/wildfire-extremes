
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
    mutate(model = burn_area_fits[i])
  train_ba_loglik[[i]] <- post$loglik_f %>%
    reshape2::melt(varnames = c('iter', 'idx')) %>%
    as_tibble %>%
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

ba_ll_df <- holdout_ba_loglik %>%
  group_by(iter, model, train) %>%
  summarize(mean_ll = mean(value))

train_ba_ll_df <- train_ba_loglik %>%
  group_by(iter, model, train) %>%
  summarize(mean_ll = mean(value))

# data frame for plotting colors
cols <- c('Gamma' = 'purple',
          'Weibull' = 'darkgreen',
          'Generalized Pareto' = 'red',
          'Tapered Pareto' = 'orange',
          'Lognormal' = 'dodgerblue')


ba_ll_df %>%
  full_join(train_ba_ll_df) %>%
  ungroup %>%
  mutate(train = ifelse(train == TRUE, 'train', 'test'),
         Distribution = case_when(
           grepl('gamma', .$model) ~ 'Gamma',
           grepl('lognormal', .$model) ~ 'Lognormal',
           grepl('_pareto', .$model) ~ 'Generalized Pareto',
           grepl('_tpareto', .$model) ~ 'Tapered Pareto',
           grepl('weibull', .$model) ~ 'Weibull'
         )) %>%
  spread(train, mean_ll) %>%
  ggplot(aes(x = train, y = test, color = Distribution)) +
  geom_point(alpha = .05) +
  xlab('Log likelihood: training set') +
  ylab('Log likelihood: test set') +
  scale_color_manual('Burn area distribution', values = cols) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0)) +
  ylim(c(-11.4, -9.6)) +
  theme_minimal()

## Generate plots to evaluate distributional assumptions
train_ba_rep <- train_ba_rep %>%
  bind_rows

# Whole distribution
train_ba_rep %>%
  filter(iter < 1001) %>%
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
  facet_wrap(~Distribution) +
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

train_ba_rep %>%
  filter(!grepl('gamma', model),
         !grepl('weibull', model)) %>%
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
  scale_color_manual('Burn area distribution', values = cols) +
  geom_point() +
  geom_errorbarh(aes(xmin = lo, xmax = hi)) +
  geom_errorbarh(aes(xmin = q25, xmax = q75), size = 1.25) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~Distribution) +
  coord_cartesian(xlim = c(1e5, 1e8),
                  ylim = c(1e-4, .01)) +
  theme(legend.position = 'none') +
  geom_point(data = obs_ccdf, color = 'black', aes(x = values), shape = 1) +
  xlab('Burn area exceedance (acres)') +
  ylab('Complementary cumulative distribution function')
