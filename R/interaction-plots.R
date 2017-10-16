library(raster)
library(rstan)
library(sf)
library(rmapshaper)
library(tidyverse)
library(magick)
library(snowfall)
library(loo)

source('R/02-explore.R')
m_fit <- read_rds('m_fit.rds')
post <- rstan::extract(m_fit)
rm(m_fit)
gc()

beta_df <- post$beta %>%
  reshape2::melt(varnames = c('iter', 'dim', 'col')) %>%
  tbl_df %>%
  group_by(dim, col) %>%
  summarize(median = median(value),
            lo = quantile(value, 0.05),
            hi = quantile(value, 0.95),
            p_neg = mean(value < 0),
            p_pos = mean(value > 0)) %>%
  ungroup %>%
  mutate(variable = colnamesX[col],
         nonzero = p_neg > .9 | p_pos > .9)

nz_beta <- beta_df %>%
  filter(nonzero)


# plot the contribution of each explanatory variable on a response
st_covs$row <- 1:nrow(st_covs)

mu_df <- post$mu_full %>%
  reshape2::melt(varnames = c('iter', 'response', 'row')) %>%
  tbl_df %>%
  group_by(response, row) %>%
  summarize(median = median(value),
            lo = quantile(value, 0.05),
            hi = quantile(value, 0.95)) %>%
  ungroup

## Prototype the plot on a subset of the data -------------------------------
which_dim <- 3
which_ecoregion <- sample(unique(st_covs$NA_L3NAME), 1)

subd <- mu_df %>%
  left_join(distinct(st_covs, NA_L3NAME, ym, row, year)) %>%
  filter(NA_L3NAME == which_ecoregion, response == which_dim, year < cutoff_year)

subd %>%
  ggplot(aes(ym, y = exp(median))) +
  geom_ribbon(aes(ymin = exp(lo), ymax = exp(hi)), alpha = .2) +
  geom_line() +
  xlab('Time') +
  ylab('Expected number of fires > 1000 acres') +
  scale_y_log10()

X_sub <- X[subd$row, ]

# for each variable and each ym, compute the product of the covariate and the coefficient
n_iter <- length(post$lp__)
elementwise_effs <- array(0, dim = c(nrow(X_sub), ncol(X_sub)))
colnames(elementwise_effs) <- colnames(X_sub)
assert_that(ncol(X_sub) == dim(post$beta)[3])
pb <- txtProgressBar(max = nrow(X_sub), style = 3)
for (i in 1:nrow(X_sub)) {
  for (j in 1:ncol(X_sub)) {
    if (X_sub[i, j] != 0) {
      elementwise_effs[i, j] <- median(X_sub[i, j] * post$beta[, which_dim, j])
    }
  }
  setTxtProgressBar(pb, i)
}

# subset this matrix to columns that are not all equal to zero
all_zero <- apply(elementwise_effs, 2, FUN = function(x) sum(x) == 0)
elementwise_effs <- elementwise_effs[, !all_zero]

# subset further to columns with nonzero coefficients
effects_to_plot <- elementwise_effs[, colnames(elementwise_effs) %in% nz_beta$variable] %>%
  as_tibble

bind_cols(subd, effects_to_plot) %>%
  gather(explanatory_variable, effect, -response, -row, -median,
         -lo, -hi, -NA_L3NAME, -year, -ym) %>%
  ggplot(aes(ym, y = effect)) +
  geom_line() +
  facet_wrap(~explanatory_variable) +
  xlab('Time') +
  ylab('Contribution to the expected number of fires') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  ggtitle(paste('Esitmated effects for ', which_ecoregion))
