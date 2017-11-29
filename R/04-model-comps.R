
# Create tables for LOO-IC for each fit -----------------------------------
library(loo)
library(tidyverse)
library(cowplot)
library(extraDistr)
library(ggridges)
library(ggExtra)

source('R/02-explore.R')
source('R/make-stan-d.R')

# possibly fetch from Amazon
#system("aws s3 cp s3://earthlab-mjoseph .. --recursive --exclude '*' --include '*.rds'")

model_fits <- list.files(pattern = '*.fit.*\\.rds')

burn_area_fits <- grep(model_fits, pattern = 'ba_', value = TRUE)

# for each model fit, produce a vector of the holdout log likelihoods
holdout_ba_loglik <- list()
train_ba_loglik <- list()
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
}

holdout_ba_loglik <- bind_rows(holdout_ba_loglik) %>%
  mutate(train = FALSE)

train_ba_loglik <- bind_rows(train_ba_loglik) %>%
  mutate(train = TRUE)

ba_ll_df <- holdout_ba_loglik %>%
  group_by(iter, model, train) %>%
  summarize(mean_nll = mean(-value))

train_ba_ll_df <- train_ba_loglik %>%
  group_by(iter, model, train) %>%
  summarize(mean_nll = mean(-value))

ba_ll_df %>%
  full_join(train_ba_ll_df) %>%
  ungroup %>%
  mutate(train = ifelse(train == TRUE, 'train', 'test'),
         Distribution = case_when(
           grepl('gamma', .$model) ~ 'Gamma',
           grepl('lognormal', .$model) ~ 'Lognormal',
           grepl('pareto', .$model) ~ 'Generalized Pareto',
           grepl('weibull', .$model) ~ 'Weibull'
         )) %>%
  spread(train, mean_nll) %>%
  ggplot(aes(x = train, y = test, color = Distribution)) +
  geom_point(alpha = .05) +
  xlab('Negative log likelihood: training set') +
  ylab('Negative log likelihood: test set') +
  scale_color_gdocs('Distribution for burn area') +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.justification = c(0, 1), legend.position = c(0, 1))
