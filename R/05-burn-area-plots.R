source('R/02-explore.R')
library(extraDistr)
library(ggridges)
library(ggthemes)

fit <- read_rds(path = list.files(pattern = 'lognormal_.*'))

# Evaluate convergence ----------------------------------------------------
traceplot(fit, inc_warmup = TRUE)

traceplot(fit, pars = c('tau', 'c', 'alpha'))
traceplot(fit, pars = c('Rho_beta'))

# Extract posterior draws and visualize results ---------------------------
post <- rstan::extract(fit)
str(post)
rm(fit)
gc()



# Coefficients that seem to be far from zero ------------------------------
beta_df <- post$beta %>%
  reshape2::melt(varnames = c('iter', 'dim', 'col')) %>%
  tbl_df

beta_summary <- beta_df %>%
  group_by(dim, col) %>%
  summarize(median = median(value),
            lo = quantile(value, 0.025),
            hi = quantile(value, 0.975),
            p_neg = mean(value < 0),
            p_pos = mean(value > 0)) %>%
  ungroup %>%
  mutate(variable = colnamesX[col],
         nonzero = p_neg > .8 | p_pos > .8)

beta_summary %>%
  filter(nonzero) %>%
  select(dim, col, p_neg, p_pos, variable) %>%
  left_join(beta_df) %>%
  ggplot(aes(value, variable)) +
  geom_density_ridges() +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  facet_wrap(~ dim, scales = 'free') +
  theme(axis.text.y = element_text(size = 7))

rm(beta_df)
gc()

# Visualize some predictions ----------------------------------------------
st_covs$row <- 1:nrow(st_covs)

mu_df <- post$mu_full %>%
  reshape2::melt(varnames = c('iter', 'row')) %>%
  tbl_df %>%
  group_by(row) %>%
  summarize(median = median(value),
            lo = quantile(value, 0.025),
            hi = quantile(value, 0.975),
            q1 = quantile(value, .1),
            q2 = quantile(value, .2),
            q3 = quantile(value, .3),
            q4 = quantile(value, .4),
            q6 = quantile(value, .6),
            q7 = quantile(value, .7),
            q8 = quantile(value, .8),
            q9 = quantile(value, .9)) %>%
  ungroup

## Visualize parameter time series
plot_mu_ts <- function(df) {
  df %>%
    ggplot(aes(ym, exp(median), fill = NA_L1NAME)) +
    theme_minimal() +
    geom_ribbon(aes(ymin = exp(lo), ymax = exp(hi)),
                alpha = .3, color = NA) +
    geom_ribbon(aes(ymin = exp(q1), ymax = exp(q9)),
                alpha = .3, color = NA) +
    geom_ribbon(aes(ymin = exp(q2), ymax = exp(q8)),
                alpha = .3, color = NA) +
    geom_ribbon(aes(ymin = exp(q3), ymax = exp(q7)),
                alpha = .3, color = NA) +
    geom_ribbon(aes(ymin = exp(q4), ymax = exp(q6)),
                alpha = .3, color = NA) +
    geom_line(size = .1) +
    facet_wrap(~ facet_factor) +
    xlab('Date') +
    scale_y_log10() +
    geom_vline(xintercept = cutoff_year,
               linetype = 'dashed', col = 'grey') +
    scale_color_gdocs() +
    scale_fill_gdocs() +
    theme(legend.position = 'none')
}

# location parameter
mu_df %>%
  full_join(st_covs) %>%
  left_join(dplyr::distinct(tbl_df(ecoregions),
                            NA_L3NAME, NA_L1CODE)) %>%
  mutate(facet_factor = paste(NA_L1CODE, NA_L3NAME)) %>%
  plot_mu_ts +
  ylab('Lognormal mean')


