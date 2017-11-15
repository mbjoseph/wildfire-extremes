source('R/02-explore.R')
library(extraDistr)
library(ggridges)

fit <- read_rds(path = list.files(pattern = 'zinbgpdfit_.*')[3])

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

wide_beta <- beta_summary %>%
  select(col, dim, median, variable) %>%
  spread(dim, median) %>%
  rename(Lomax_shape = `1`,
         Lomax_scale = `2`,
         NegBinom_precision = `3`,
         NegBinom_mean = `4`,
         logit_p_zero = `5`)

wide_beta %>%
  select(-col, -variable) %>%
  pairs

rm(wide_beta)
rm(beta_df)
gc()

# Visualize some predictions ----------------------------------------------
st_covs$row <- 1:nrow(st_covs)

mu_df <- post$mu_full %>%
  reshape2::melt(varnames = c('iter', 'response', 'row')) %>%
  tbl_df %>%
  group_by(response, row) %>%
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

# shape parameter
mu_df %>%
  filter(response == 1) %>%
  full_join(st_covs) %>%
  left_join(dplyr::distinct(tbl_df(ecoregions),
                            NA_L3NAME, NA_L1CODE)) %>%
  mutate(facet_factor = paste(NA_L1CODE, NA_L3NAME)) %>%
  plot_mu_ts +
  ylab('Pareto 2 shape parameter')

# scale parameter
mu_df %>%
  filter(response == 2) %>%
  full_join(st_covs) %>%
  left_join(dplyr::distinct(tbl_df(ecoregions),
                            NA_L3NAME, NA_L1CODE)) %>%
  mutate(facet_factor = paste(NA_L1CODE, NA_L3NAME)) %>%
  plot_mu_ts +
  ylab('Pareto 2 scale parameter')

# neg. binom precision
mu_df %>%
  filter(response == 3) %>%
  full_join(st_covs) %>%
  left_join(dplyr::distinct(tbl_df(ecoregions),
                            NA_L3NAME, NA_L1CODE)) %>%
  mutate(facet_factor = paste(NA_L1CODE, NA_L3NAME)) %>%
  plot_mu_ts +
  ylab('Negative binomial precision')

# neg. binom mean
mu_df %>%
  filter(response == 4) %>%
  full_join(st_covs) %>%
  left_join(dplyr::distinct(tbl_df(ecoregions),
                            NA_L3NAME, NA_L1CODE)) %>%
  mutate(facet_factor = paste(NA_L1CODE, NA_L3NAME)) %>%
  plot_mu_ts +
  ylab('Negative binomial mean')

mu_df %>%
  filter(response == 5) %>%
  full_join(st_covs) %>%
  left_join(dplyr::distinct(tbl_df(ecoregions),
                            NA_L3NAME, NA_L1CODE)) %>%
  mutate(facet_factor = paste(NA_L1CODE, NA_L3NAME)) %>%
  ggplot(aes(ym, plogis(median), fill = NA_L1NAME)) +
  theme_minimal() +
  geom_ribbon(aes(ymin = plogis(lo), ymax = plogis(hi)),
              alpha = .3, color = NA) +
  geom_ribbon(aes(ymin = plogis(q1), ymax = plogis(q9)),
              alpha = .3, color = NA) +
  geom_ribbon(aes(ymin = plogis(q2), ymax = plogis(q8)),
              alpha = .3, color = NA) +
  geom_ribbon(aes(ymin = plogis(q3), ymax = plogis(q7)),
              alpha = .3, color = NA) +
  geom_ribbon(aes(ymin = plogis(q4), ymax = plogis(q6)),
              alpha = .3, color = NA) +
  geom_line(size = .1) +
  facet_wrap(~ facet_factor) +
  xlab('Date') +
  geom_vline(xintercept = cutoff_year,
             linetype = 'dashed', col = 'grey') +
  scale_color_gdocs() +
  scale_fill_gdocs() +
  theme(legend.position = 'none') +
  ylab('Probability of excess zero')
ggsave(filename = 'fig/p-excess-zero.pdf', width = 30, height = 10)

# Compare predicted to expected counts for training data
mu_df %>%
  filter(response == 4) %>%
  full_join(st_covs) %>%
  left_join(dplyr::distinct(tbl_df(ecoregions),
                            NA_L3NAME, NA_L1CODE)) %>%
  mutate(facet_factor = paste(NA_L1CODE, NA_L3NAME)) %>%
  dplyr::select(ym, NA_L3NAME, median, lo, hi) %>%
  right_join(count_df) %>%
  full_join(er_df) %>%
  mutate(is_train = year < cutoff_year) %>%
  filter(!is_train, n_fire > 0) %>%
  ggplot(aes(x = n_fire, y = exp(median), color = is_train)) +
  facet_wrap(~ NA_L1NAME) +
  geom_point(alpha = .5) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'blue') +
  geom_segment(aes(xend = n_fire, y = exp(lo), yend = exp(hi)), alpha = .5) +
  scale_x_log10() +
  scale_y_log10(limits = c(.001, 1e4))

gc()

