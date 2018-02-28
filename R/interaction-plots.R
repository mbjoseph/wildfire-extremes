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
library(ggthemes)

source('R/02-explore.R')
m_fit <- read_rds('zinb_full_fit.rds')
post <- rstan::extract(m_fit)
ba_post <- rstan::extract(read_rds('ba_lognormal_fit.rds'))
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

mtbs_totals <- mtbs %>%
  select(-geometry) %>%
  as_tibble %>%
  group_by(ym, NA_L3NAME) %>%
  summarize(total_area = sum(R_ACRES))

unique_ecoregions <- unique(st_covs$NA_L3NAME)


write_attribution_plot <- function(which_ecoregion,
                                   max_iter = 1000) {

  which_ecoregion <- 'Arizona/New Mexico Mountains'
  max_iter <- 100

  plot_name <- file.path('fig', 'effs',
                         paste0(tolower(gsub(' |/', '-', which_ecoregion)),
                                '.png'))

  subd <- mu_df %>%
    left_join(distinct(st_covs, NA_L3NAME, ym, row, year)) %>%
    filter(NA_L3NAME == which_ecoregion,
           year >= cutoff_year,
           response == 1) %>%
    mutate(row_index = 1:n())

  # subset design matrix
  X_sub <- X[subd$row, ]
  nonzero_columns <- which(apply(X_sub, 2, function(x) !all(x == 0)))
  X_sub <- X_sub[, nonzero_columns]

  # for each response, variable, and ym, multiply covariate and coefficient
  n_iter <- min(max_iter, length(post$lp__))
  elementwise_effs <- array(0, dim = c(nrow(X_sub), ncol(X_sub), n_iter, 2))
  colnames(elementwise_effs) <- colnames(X_sub)
  effects_to_plot <- list()
  pb <- txtProgressBar(max = nrow(X_sub), style = 3)
  for (i in 1:nrow(X_sub)) {
    for (j in 1:ncol(X_sub)) {
      for (response in 1:2) {
        if (X_sub[i, j] != 0) {
          elementwise_effs[i, j, , response] <- X_sub[i, j] *
                      post$beta[1:max_iter, response, nonzero_columns[j]]
        }
      }
    }
    effects_to_plot[[i]] <- elementwise_effs[i, , , ] %>%
      reshape2::melt(varnames = c('col', 'iter', 'response')) %>%
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
             row_index = i) %>%
      group_by(row_index, variable, iter, response) %>%
      summarize(total_effect = sum(value))
    setTxtProgressBar(pb, i)
  }

  # aggregate by variable
  effects_summary <- effects_to_plot %>%
    bind_rows %>%
    group_by(row_index, variable, response) %>%
    summarize(median_eff = median(total_effect),
              lo_eff = quantile(total_effect, .05),
              hi_eff = quantile(total_effect, .95)) %>%
    ungroup %>%
    left_join(select(subd, row_index, ym, NA_L3NAME, year)) %>%
    ungroup %>%
    mutate(variable = tolower(variable),
           variable = tools::toTitleCase(variable))


  # median plot
  effects_summary %>%
    group_by(variable, response) %>%
    filter(!grepl('^L', variable),
           variable != 'Terrain Ruggedness',
           year >= cutoff_year) %>%
    ggplot(aes(ym, y = median_eff, color = variable)) +
    geom_line() +
    geom_point(size = .2) +
    facet_wrap(~response, scales = 'free_y', ncol = 1) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank()) +
    xlab('Time') +
    ylab('Contribution to fire risk') +
    ggtitle(which_ecoregion) +
    scale_color_discrete('Variable')
  med_plot_name <- file.path('fig', 'effs',
                         paste0(tolower(gsub(' |/', '-', which_ecoregion)),
                                '-', which_dim, '-', 'med',
                                '.png'))
  ggsave(filename = med_plot_name, plot = p_med, width = 8, height = 3.5)
  return(list(plot = p_med, effects = effects_to_plot))
}

# choose a case study: why not the biggest fire in the train set?
max_d <- holdout_burns[which.max(holdout_burns$R_ACRES), ]

write_attribution_plot(max_d$NA_L3NAME, which_dim = 1)


# Create for every ecoregion (takes a long time)
# dir.create(file.path('fig', 'effs'))
#
# plot_list <- list()
# for (i in seq_along(unique_ecoregions)) {
#   plot_list[[i]] <- write_attribution_plot(unique_ecoregions[i], which_dim = 1)
# }

