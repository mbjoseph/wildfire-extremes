
# Model comparisons for counts --------------------------------------------
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

count_fits <- grep(model_fits, pattern = 'ba_', value = TRUE, invert = TRUE)

# for each model fit, produce a vector of the holdout log likelihoods
holdout_c_loglik <- list()
#train_c_loglik <- list()
for (i in seq_along(count_fits)) {
  post <- rstan::extract(read_rds(count_fits[i]))
  holdout_c_loglik[[i]] <- post$holdout_loglik_c %>%
    reshape2::melt(varnames = c('iter', 'idx')) %>%
    as_tibble %>%
    mutate(model = count_fits[i])
  # train_c_loglik[[i]] <- post$train_loglik_c %>%
  #   reshape2::melt(varnames = c('iter', 'idx')) %>%
  #   as_tibble %>%
  #   mutate(model = burn_area_fits[i])
}

holdout_c_loglik <- bind_rows(holdout_c_loglik) %>%
  mutate(train = FALSE)

# train_c_loglik <- bind_rows(train_c_loglik) %>%
#   mutate(train = TRUE)

c_ll_df <- holdout_c_loglik %>%
  group_by(iter, model, train) %>%
  summarize(mean_nll = mean(-value))

# train_c_ll_df <- train_c_loglik %>%
#   group_by(iter, model, train) %>%
#   summarize(mean_nll = mean(-value))

c_ll_df %>%
  ungroup %>%
  mutate(Distribution = case_when(
    grepl('nb', .$model) ~ 'Negative binomial',
    grepl('pois', .$model) ~ 'Poisson'
  )) %>%
  ggplot(aes(x = -mean_nll, fill = Distribution)) +
  geom_density(alpha = .5) +
  scale_fill_gdocs('Distribution for counts') +
  xlab('Log likelihood: test set') +
  ylab('Posterior density') +
  xlim(-.9, -.44) +
  theme(legend.justification = c(0, 1), legend.position = c(0, 1))
# TODO: re-run models with training loglik calculations and make scatterplot
