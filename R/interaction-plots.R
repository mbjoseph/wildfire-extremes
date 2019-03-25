library(raster)
library(rstan)
library(sf)
library(tidyverse)
library(ggthemes)
library(patchwork)
library(assertthat)

st_covs <- read_rds('data/processed/st_covs.rds')
holdout_burns <- read_rds('data/processed/holdout_burns.rds')
cutoff_year <- read_rds('data/processed/cutoff_year.rds')
X <- read_rds('data/processed/X.rds')
mtbs <- read_rds('data/processed/mtbs.rds')
stan_d <- read_rds('data/processed/stan_d.rds')
min_size <- stan_d$min_size

post <- rstan::extract(read_rds('zinb_full_fit.rds'), 
                       pars = c('mu_full', 'lp__', 'beta'))
ba_post <- rstan::extract(read_rds('ba_lognormal_fit.rds'), 
                          pars = 'beta')
gc()

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


write_attribution_plot <- function(which_ecoregion) {
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
  n_iter <- length(post$lp__)
  # 3 dims, one for each response (2 for count, one for burn area)
  elementwise_effs <- array(0, dim = c(nrow(X_sub), ncol(X_sub), n_iter, 3))
  colnames(elementwise_effs) <- colnames(X_sub)
  effects_to_plot <- list()
  pb <- txtProgressBar(max = nrow(X_sub), style = 3)
  for (i in 1:nrow(X_sub)) {
    for (j in 1:ncol(X_sub)) {
      for (response in 1:3) {
        if (X_sub[i, j] != 0) {
          if (response == 1) {
            # burn area mean
            elementwise_effs[i, j, , response] <- X_sub[i, j] *
              ba_post$beta[, response, nonzero_columns[j]]
          } else {
            # count params
            elementwise_effs[i, j, , response] <- X_sub[i, j] *
              post$beta[, response - 1, nonzero_columns[j]]
          }
        }
      }
    }
    effects_to_plot[[i]] <- elementwise_effs[i, , , ] %>%
      reshape2::melt(varnames = c('col', 'iter', 'response')) %>%
      as_tibble %>%
      mutate(col = as.character(col),
             variable = case_when(grepl('tri', .$col) ~ 'Terrain ruggedness',
                                  grepl('log_housing_density', .$col) ~ 'Housing density',
                                  grepl('vs', .$col) ~ 'Wind speed',
                                  grepl('pr_', .$col) ~ 'Monthly precip.',
                                  grepl('prev_12mo_precip', .$col) ~ '12 month precip.',
                                  grepl('tmmx', .$col) ~ 'Temperature',
                                  grepl('rmin', .$col) ~ 'Humidity'),
             variable = ifelse(is.na(variable), col, variable),
             variable = gsub('NA_', '', variable),
             variable = gsub('NAME', ' ', variable),
             variable = gsub('_', ': ', variable),
             row_index = i) %>%
      group_by(row_index, variable, iter, response) %>%
      summarize(total_effect = sum(value))
    assert_that(!any(is.na(effects_to_plot[[i]]$variable)))
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
  plot_data <- effects_summary %>%
    filter(!grepl('^L', variable),
           variable != 'Terrain Ruggedness',
           year >= cutoff_year) %>%
    mutate(Response = case_when(.$response == 1 ~ 'Lognormal mean (burn area)',
                                .$response == 2 ~ 'Negative binomial mean (number of fires)',
                                .$response == 3 ~ 'Zero-inflation component (number of fires)'),
           variable = factor(tools::toTitleCase(tolower(variable)),
                             levels = c('Humidity', 'Temperature',
                                        'Housing Density', 'Monthly Precip.',
                                        '12 Month Precip.', 'Wind Speed')))

  p_med <- plot_data %>%
    ggplot(aes(ym, y = median_eff, color = variable)) +
    geom_line() +
    facet_wrap(~Response, scales = 'free_y', ncol = 1) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank()) +
    xlab('') +
    ylab('Contribution to linear predictor') +
    scale_color_manual('', values = c('dodgerblue',
                                              'red',
                                              scales::alpha(
                                                c('green4',
                                                  'pink',
                                                  'orange',
                                                  'lightblue'), .7)))
  return(list(plot = p_med,
              effects = effects_to_plot %>%
                bind_rows %>%
                left_join(select(subd, row_index, ym, NA_L3NAME, year)),
              plot_data = plot_data))
}

# choose a case study: why not the biggest fire in the train set?
max_d <- holdout_burns[which.max(holdout_burns$Acres), ]

case_study_plot <- write_attribution_plot(max_d$NA_L3NAME)

p1 <- case_study_plot$plot +
  geom_vline(aes(xintercept = ym), data = max_d, lty = 'dotted')
p1
ggsave(filename = 'fig/figure_11.pdf', plot = p1,
       width = 7, height = 4)
write_rds(case_study_plot, 'data/processed/case-study-plot.rds')
