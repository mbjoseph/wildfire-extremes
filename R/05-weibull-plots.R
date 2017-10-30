source('R/02-explore.R')

w_fit <- read_rds(path = list.files(pattern = 'wfit_.*'))

# Evaluate convergence ----------------------------------------------------
traceplot(w_fit, inc_warmup = TRUE)

traceplot(w_fit, pars = c('tau', 'c', 'alpha'))
traceplot(w_fit, pars = c('Rho_beta'))

plot(w_fit, pars = 'beta') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  coord_flip()


# Extract posterior draws and visualize results ---------------------------
post <- rstan::extract(w_fit)
str(post)
rm(w_fit)
gc()



# Coefficients that seem to be far from zero ------------------------------
beta_df <- post$beta %>%
  reshape2::melt(varnames = c('iter', 'dim', 'col')) %>%
  tbl_df %>%
  group_by(dim, col) %>%
  summarize(median = median(value),
            lo = quantile(value, 0.025),
            hi = quantile(value, 0.975),
            p_neg = mean(value < 0),
            p_pos = mean(value > 0)) %>%
  ungroup %>%
  mutate(variable = colnamesX[col],
         nonzero = p_neg > .95 | p_pos > .95)

beta_df %>%
  filter(nonzero) %>%
  ggplot(aes(median, variable)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_segment(aes(x = lo, xend = hi, yend = variable)) +
  facet_wrap(~ dim, scales = 'free') +
  theme(axis.text.y = element_text(size = 7))

wide_beta <- beta_df %>%
  select(col, dim, median, variable) %>%
  spread(dim, median) %>%
  rename(Weibull_shape = `1`,
         Weibull_scale = `2`,
         NegBinom_precision = `3`,
         NegBinom_mean = `4`)

wide_beta %>%
  select(-col, -variable) %>%
  pairs


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
  ylab('Weibull shape parameter')

# scale parameter
mu_df %>%
  filter(response == 2) %>%
  full_join(st_covs) %>%
  left_join(dplyr::distinct(tbl_df(ecoregions),
                            NA_L3NAME, NA_L1CODE)) %>%
  mutate(facet_factor = paste(NA_L1CODE, NA_L3NAME)) %>%
  plot_mu_ts +
  ylab('Weibull scale parameter')

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
  filter(!is_train) %>%
  ggplot(aes(x = n_fire, y = exp(median), color = is_train)) +
  facet_wrap(~ NA_L1NAME) +
  geom_point(alpha = .5) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'blue') +
  geom_segment(aes(xend = n_fire, y = exp(lo), yend = exp(hi)), alpha = .5) +
  scale_x_log10() +
  scale_y_log10(limits = c(.001, 1e4))

gc()

## Posterior predictive checks ---------------------------

# set breakpoints for histogram
max_breaks <- 200
mtbs_hist <- hist(mtbs$R_ACRES, breaks = max_breaks, right = FALSE)
write_rds(mtbs_hist, 'mtbs_hist.rds')

n_draw <- dim(post$mu_full)[1]
n_preds <- dim(post$mu_full)[3]
hist_counts <- matrix(nrow = n_draw, ncol = max_breaks*1.5)
hist_mids <- matrix(nrow = n_draw, ncol = max_breaks*1.5)
block_maxima <- matrix(nrow = n_draw, ncol = n_preds)
area_sums <- matrix(nrow = n_draw, ncol = n_preds)
n_events <- rnbinom(n_draw * n_preds,
                    size = exp(c(post$mu_full[, 3, ])),
                    mu = exp(c(post$mu_full[, 4, ]))) %>%
  matrix(nrow = n_draw, ncol = n_preds)
sum(is.na(n_events))
gc()

pb <- txtProgressBar(max = n_draw, style = 3)
for (i in 1:n_draw) {
  ba_vec <- list() # to hold all fire sizes for iteration i
  counter <- 1
  for (j in 1:n_preds) {
    if (!is.na(n_events[i, j]) & n_events[i, j] > 0) {
      burn_areas <- 1e3 +
        weibull_scale_adj * rweibull(n = n_events[i, j],
                                   shape = exp(post$mu[i, 1, j]),
                                   scale = exp(post$mu[i, 2, j]))
      block_maxima[i, j] <- max(burn_areas)
      area_sums[i, j] <- sum(burn_areas)
      ba_vec[[counter]] <- burn_areas
      counter <- counter + 1
    }
  }
  ba_vec <- unlist(ba_vec)
  ba_hist <- hist(ba_vec, breaks = max_breaks, plot = FALSE, right = FALSE)
  nbins <- length(ba_hist$counts)
  hist_counts[i, 1:nbins] <- ba_hist$counts
  hist_mids[i, 1:nbins] <- ba_hist$mids
  setTxtProgressBar(pb, i)
}

ppred_df <- reshape2::melt(block_maxima,
                           varnames = c('iter', 'idx'),
                           value.name = 'max_burn_area') %>%
  bind_cols(reshape2::melt(area_sums,
                           varnames = c('iter', 'idx'),
                           value.name = 'total_burn_area')) %>%
  bind_cols(reshape2::melt(n_events,
                           varnames = c('iter', 'idx'),
                           value.name = 'n_events')) %>%
  as_tibble %>%
  dplyr::select(-ends_with('1'), -ends_with('2')) %>%
  filter(!is.na(n_events))

write_rds(ppred_df, 'ppred_df.rds')
rm(ppred_df)
gc()

hist_df <- reshape2::melt(hist_counts,
                          varnames = c('iter', 'bin'),
                          value.name = 'count') %>%
  bind_cols(reshape2::melt(hist_mids,
                           varnames = c('iter', 'bin'),
                           value.name = 'midpoint')) %>%
  as_tibble %>%
  select(-ends_with('1')) %>%
  arrange(iter, bin) %>%
  na.omit %>%
  mutate(approx_burn_area = count * midpoint) %>%
  group_by(iter) %>%
  mutate(cum_burn_area = cumsum(approx_burn_area)) %>%
  ungroup

write_rds(hist_df, 'ppc_hist.rds')

