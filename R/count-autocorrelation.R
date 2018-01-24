
# Check for autocorrelation in residuals --------------------------
library(raster)
library(tidyverse)
library(ggridges)
library(patchwork)
library(viridis)
library(purrr)

source('R/02-explore.R')
source('R/make-stan-d.R')

zinb_post <- rstan::extract(read_rds('zinb_fit.rds'))

# compute pearson residuals
assert_that(all(diff(match(count_df$er_ym, st_covs$er_ym)) == 1))

pearson_resids <- matrix(nrow = nrow(count_df), ncol = length(zinb_post$lp__))
pb <- txtProgressBar(max = ncol(pearson_resids))
for (i in 1:ncol(pearson_resids)) {
  pearson_resids[, i] <- (count_df$n_fire - zinb_post$zinb_mean[i, ]) /
    sqrt(zinb_post$zinb_var[i, ])
  setTxtProgressBar(pb, i)
}

resid_df <- count_df %>%
  dplyr::select(NA_L3NAME, ym, year, er_ym) %>%
  mutate(idx = 1:n())

resid_df <- pearson_resids %>%
  reshape2::melt(varnames = c('idx', 'iter'), value.name = 'pearson_resid') %>%
  as_tibble %>%
  left_join(resid_df) %>%
  filter(year < 2010) %>%
  mutate(er_iter = paste(NA_L3NAME, iter, sep = '_'))

resid_df

acf_df <- resid_df %>%
  filter(iter < 100) %>%
  dplyr::select(er_iter, pearson_resid) %>%
  split(.$er_iter) %>%
  map(~acf(.$pearson_resid, lag.max = 2, plot = FALSE)) %>%
  map(~tibble(acf = .$acf[-1], lag = .$lag[-1])) %>%
  bind_rows(.id = 'er_iter') %>%
  separate(er_iter, into = c('NA_L3NAME', 'iter'), sep = '_')

acf_df %>%
  filter(lag == 1) %>%
  group_by(NA_L3NAME) %>%
  mutate(p_pos = mean(acf > 0),
         med_acf = median(acf)) %>%
  ggplot(aes(x = acf, y = reorder(NA_L3NAME, med_acf), fill = p_pos)) +
  geom_density_ridges() +
  scale_fill_viridis_c()


# visualize time series for random effects
phi_df <- zinb_post$phi %>%
  reshape2::melt(varnames = c('iter', 'm', 'n', 't')) %>%
  group_by(m, n, t) %>%
  summarize(med = median(value),
            lo = quantile(value, .1),
            hi = quantile(value, .9))

phi_df %>%
  filter(m == 2) %>%
  ggplot(aes(t, med, fill = factor(m), color = factor(m))) +
  geom_ribbon(aes(ymin = lo, ymax = hi), color = NA, alpha = .5) +
  geom_line() +
  facet_wrap(~n, scales = 'free_y')
