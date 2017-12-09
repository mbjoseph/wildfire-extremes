
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

count_fits <- grep(model_fits, pattern = 'ba_', value = TRUE, invert = TRUE)

# for each model fit, produce a vector of the holdout log likelihoods
holdout_c_loglik <- list()
train_c_loglik <- list()
train_c_rep <- list()
holdout_c_rep <- list()

for (i in seq_along(count_fits)) {
  post <- rstan::extract(read_rds(count_fits[i]))
  holdout_c_loglik[[i]] <- post$holdout_loglik_c %>%
    apply(1, sum) %>%
    reshape2::melt(varnames = 'iter', value.name = 'holdout_loglik')
}

names(holdout_c_loglik) <- count_fits

holdout_c_loglik <- bind_rows(holdout_c_loglik, .id = 'model') %>%
  as_tibble()

holdout_c_loglik %>%
  ggplot(aes(x = holdout_loglik, y = model)) +
  geom_jitter(alpha = .1)
