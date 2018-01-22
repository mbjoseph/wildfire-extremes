library(raster)
library(rstan)
library(sf)
library(rmapshaper)
library(tidyverse)
library(magick)
library(snowfall)
library(loo)
library(ggExtra)

source('R/02-explore.R')
m_fit <- read_rds('zinb_fit.rds')
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
         nonzero = p_neg > .95 | p_pos > .95)

nz_beta <- beta_df %>%
  filter(nonzero)

med_df <- beta_df %>%
  select(dim, median, variable) %>%
  spread(dim, median) %>%
  rename(nb_eff = `1`,
         p_eff = `2`)

med_df %>%
  ggplot(aes(nb_eff, p_eff)) +
  geom_hline(yintercept = 0, alpha = .2) +
  geom_vline(xintercept = 0, alpha = .2) +
  geom_point() +
  geom_bin2d(binwidth = c(.004, .004), color = NA) +
  theme_minimal() +
  scale_fill_viridis_c(trans = 'log', option = 'B',
                       breaks = c(1, 10, 100, 1000), 'Count') +
  xlab('Coefficient: negative binomial component') +
  ylab('Coefficient: zero-inflation component')
ggsave('fig/bivariate-horseshoe.png', width = 6, height = 4)

beta_df %>%
  filter(median > .5)


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


write_attribution_plot <- function(which_ecoregion, which_dim = 1) {
  subd <- mu_df %>%
    left_join(distinct(st_covs, NA_L3NAME, ym, row, year)) %>%
    filter(NA_L3NAME == which_ecoregion, response == which_dim)

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

  p <- bind_cols(subd, effects_to_plot) %>%
    gather(explanatory_variable, effect, -response, -row, -median,
           -lo, -hi, -NA_L3NAME, -year, -ym) %>%
    ggplot(aes(ym, y = effect)) +
    geom_line() +
    facet_wrap(~explanatory_variable, ncol = 1) +
    xlab('Time') +
    ylab('Contribution to the expected number of fires') +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    ggtitle(paste('Esitmated effects for ', which_ecoregion))
  plot_name <- file.path('fig', 'effs',
                         paste0(tolower(gsub(' |/', '-', which_ecoregion)),
                                '-', which_dim,
                                '.png'))
  ggsave(filename = plot_name, plot = p, width = 10, height = 20)
}

dir.create(file.path('fig', 'effs'))

unique_ecoregions <- unique(st_covs$NA_L3NAME)
for (i in seq_along(unique_ecoregions)) {
  write_attribution_plot(unique_ecoregions[i], which_dim = 1)
}
