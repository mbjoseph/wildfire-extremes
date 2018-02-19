library(raster)
library(rstan)
library(sf)
library(rmapshaper)
library(tidyverse)
library(magick)
library(snowfall)
library(loo)
library(ggExtra)
library(viridis)

source('R/02-explore.R')
m_fit <- read_rds('zinb_full_fit.rds')
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
  geom_bin2d(binwidth = c(.1, .1), color = NA) +
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


unique_ecoregions <- unique(st_covs$NA_L3NAME)

write_attribution_plot <- function(which_ecoregion,
                                   which_dim = 1,
                                   effect_threshold = .5) {
  subd <- mu_df %>%
    left_join(distinct(st_covs, NA_L3NAME, ym, row, year)) %>%
    filter(NA_L3NAME == which_ecoregion, response == which_dim) %>%
    mutate(row_index = 1:n())

  X_sub <- X[subd$row, ]

  nonzero_columns <- which(apply(X_sub, 2, function(x) !all(x == 0)))

  X_sub <- X_sub[, nonzero_columns]

  # for each variable and each ym, compute the product of the covariate and the coefficient
  n_iter <- length(post$lp__)
  elementwise_effs <- array(0, dim = c(nrow(X_sub), ncol(X_sub)))
  colnames(elementwise_effs) <- colnames(X_sub)
  pb <- txtProgressBar(max = nrow(X_sub), style = 3)
  for (i in 1:nrow(X_sub)) {
    for (j in 1:ncol(X_sub)) {
      if (X_sub[i, j] != 0) {
        elementwise_effs[i, j] <- median(X_sub[i, j] * post$beta[, which_dim, nonzero_columns[j]])
      }
    }
    setTxtProgressBar(pb, i)
  }

  # aggregate by variable
  effects_to_plot <- elementwise_effs %>%
    reshape2::melt(varnames = c('row_index', 'col')) %>%
    as_tibble %>%
    mutate(col = as.character(col),
           variable = case_when(grepl('ctri', .$col) ~ 'Terrain ruggedness',
                                grepl('chd', .$col) ~ 'Housing density',
                                grepl('cvs', .$col) ~ 'Wind speed',
                                grepl('cpr_', .$col) ~ 'Monthly precip.',
                                grepl('cpr12_', .$col) ~ '12 month precip.',
                                grepl('ctmx', .$col) ~ 'Temperature',
                                grepl('crmin', .$col) ~ 'Humidity'),
           variable = ifelse(is.na(variable), col, variable),
           variable = gsub('NA_', '', variable),
           variable = gsub('NAME', ' ', variable),
           variable = gsub('_', ': ', variable),
           variable = tolower(variable),
           variable = tools::toTitleCase(variable)) %>%
    group_by(row_index, variable) %>%
    summarize(n = n(),
              total_effect = sum(value)) %>%
    left_join(subd) %>%
    ungroup

  p <- effects_to_plot %>%
    group_by(variable) %>%
    mutate(max_effect = max(abs(total_effect))) %>%
    filter(max_effect > effect_threshold) %>%
    ggplot(aes(ym, y = total_effect, color = variable)) +
    geom_line() +
    theme_minimal() +
    theme(panel.grid.minor = element_blank()) +
    xlab('Time') +
    ylab('Contribution to E(n)') +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    ggtitle(which_ecoregion) +
    scale_color_discrete('Variable')
  plot_name <- file.path('fig', 'effs',
                         paste0(tolower(gsub(' |/', '-', which_ecoregion)),
                                '-', which_dim,
                                '.png'))
  ggsave(filename = plot_name, plot = p, width = 6, height = 3.5)
  p
}

dir.create(file.path('fig', 'effs'))

plot_list <- list()
for (i in seq_along(unique_ecoregions)) {
  plot_list[[i]] <- write_attribution_plot(unique_ecoregions[i], which_dim = 1)
}
