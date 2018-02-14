library(tidyverse)
library(rstan)
library(viridis)
library(ggridges)
library(ggrepel)
library(hrbrthemes)
library(patchwork)
source('R/02-explore.R')
source('R/make-stan-d.R')


# Simulate from the joint predictive distribution -------------------------
# to demonstrate predictions over extremes
# using the ZINB + lognormal model

# we need the predicted counts from the zinb model
zinb_preds <- read_rds('holdout_c_rep.rds') %>%
  filter(model == 'zinb_fit.rds') %>%
  full_join(read_rds('train_c_rep.rds') %>% filter(model == 'zinb_fit.rds')) %>%
  select(-ym, -model) %>%
  left_join(select(st_covs, id, NA_L3NAME, ym)) %>%
  rename(n_event = value) %>%
  select(-year) %>%
  arrange(iter, id)

# now, having the predicted counts, we need to simulate fire sizes for
# each event, using the predictive distribution from the lognormal model
ln_fit <- read_rds('ba_lognormal_fit.rds')
ln_post <- rstan::extract(ln_fit, pars = c('mu_full', 'scale'))

# we only really need the predictions for expected value and standard
# deviation from the lognormal model, and we need to have the same
# number of posterior draws as are present in the count model
count_post_iters <- max(zinb_preds$iter)

ln_mu <- ln_post$mu_full[1:count_post_iters, ] %>%
  reshape2::melt(varnames = c('iter', 'id'), value.name = 'ln_mu') %>%
  as_tibble %>%
  arrange(iter, id)

ln_sigma <- ln_post$scale[1:count_post_iters] %>%
  reshape2::melt(varnames = c('iter'), value.name = 'ln_scale') %>%
  as_tibble

# merge lognormal mean and scale draws with counts
test_preds <- bind_cols(zinb_preds, ln_mu) %>%
  left_join(ln_sigma)

assert_that(all(test_preds$iter == test_preds$iter1))
assert_that(all(test_preds$id == test_preds$id1))

test_preds <- test_preds %>%
  select(-ends_with('1'))

# Generate derived parameters about the distribution of maxima
nmax_q <- function(p, n, mu, sigma) {
  erf.inv <- function(x) qnorm((x + 1)/2)/sqrt(2)

  sigma * sqrt(2) * erf.inv(2 * p^(1/n) - 1) + mu
}

test_preds <- test_preds %>%
  mutate(log_p = n_event * pnorm(log(100000), mean = ln_mu, sd = ln_scale, log.p = TRUE),
         q50 = nmax_q(.5, n = n_event, mu = ln_mu, sigma = ln_scale),
         q95 = nmax_q(.95, n = n_event, mu = ln_mu, sigma = ln_scale),
         q05 = nmax_q(.05, n = n_event, mu = ln_mu, sigma = ln_scale))

test_preds %>%
  filter(n_event > 0) %>%
  ggplot(aes(n_event, 1 - exp(log_p))) +
  geom_boxplot(aes(group = cut_width(n_event, 1)),
               outlier.size = .1) +
  scale_x_log10()

test_preds %>%
  filter(n_event > 0) %>%
  ggplot(aes(n_event, exp(q50) + min_size)) +
  geom_boxplot(aes(group = cut_width(n_event, 1)),
               color = 'darkred', outlier.color = NA) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank()) +
  xlab('Number of fire events') +
  ylab('50% quantile: maximum burn area (acres)') +
  xlim(0, 50) +
  scale_y_log10()

test_preds %>%
  filter(n_event < 51, n_event > 0) %>%
  ggplot(aes(exp(q95) + min_size, factor(n_event))) +
  geom_density_ridges(scale = 4, rel_min_height = 0.05,
                      fill = 'grey92', bandwidth = .12) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank()) +
  ylab('Number of fires > 1000 acres in one month') +
  xlab('95% quantile: maximum burn area (acres)') +
  scale_y_discrete(breaks = as.character(seq(0, 50, by = 10))) +
  scale_x_log10(breaks = c(1e3, 1e4, 1e5, 1e6),
                labels = c('1,000', '10,000', '100,000', '1,000,000')) +
  coord_flip()


holdout_ids <- st_covs$id[st_covs$year >= cutoff_year]
pred_df <- expand.grid(x = 1e6 - min_size,
                       id = holdout_ids,
                       iter = 1:max(test_preds$iter))

iter_max <- 600
exceedance_df <- test_preds %>%
  filter(iter < iter_max, id %in% holdout_ids) %>%
  full_join(filter(pred_df, iter < iter_max))

exceedance_summary <- exceedance_df %>%
  # find the cdf for monthly maxima
  mutate(log_p = n_event * pnorm(log(x),
                                 mean = ln_mu,
                                 sd = ln_scale,
                                 log.p = TRUE)) %>%
  left_join(select(st_covs, id, year, rmin)) %>%
  group_by(NA_L3NAME, iter, x) %>%
  # sum to get the joint probability across all time steps
  summarize(total_p = sum(log_p),
            total_events = sum(n_event),
            rmin = mean(rmin)) %>%
  group_by(NA_L3NAME, x, rmin) %>%
  summarize(`Exceedance probability` = 1 - mean(exp(total_p)),
            ep_lo = 1 - quantile(exp(total_p), .025),
            ep_hi = 1 - quantile(exp(total_p), .975),
            `E(n)` = mean(total_events),
            n_lo = quantile(total_events, .025),
            n_hi = quantile(total_events, .975),
            mean_rmin = mean(rmin)) %>%
  ungroup

point_size <- 4
n_vs_exc <- exceedance_summary %>%
  left_join(area_df) %>%
  ggplot(aes(`E(n)`, `Exceedance probability`, color = rmin)) +
  geom_point(size = point_size) +
  geom_point(shape = 1, color = 1, alpha = .5, size = point_size) +
  scale_y_log10() +
  scale_x_log10() +
  theme_minimal() +
  theme(panel.grid.minor = element_blank()) +
  xlab('Expected number of fires: 2010-2015') +
  ylab('Million acre exceedance probability') +
  scale_color_viridis(direction = -1, 'Mean humidity')
ggsave(plot = n_vs_exc,
       'fig/number-vs-exceedance.png',
       width = 7, height = 3.5)

# show on map
class(ecoregions)

ecoregions$exceedance_prob <- exceedance_summary$`Exceedance probability`[match(ecoregions$NA_L3NAME, exceedance_summary$NA_L3NAME)]

simpler_ecoregions <- ecoregions %>%
  as('Spatial') %>%
  rmapshaper::ms_simplify(keep = .01) %>%
  sf::st_as_sf()

exc_map <- simpler_ecoregions %>%
  ggplot() +
  geom_sf(aes(fill = exceedance_prob),
          color = 'white',
          size = .1) +
  scale_fill_viridis(trans = 'log', option = 'B',
                     'Million acre\nexceedance\nprobability',
                     breaks = c(.00001, .0001, .001, .01)) +
  geom_sf(data = mtbs, color = 'black', alpha = .1, size = .01) +
  hrbrthemes::theme_ipsum_rc(base_family = 'helvetica')
ggsave(plot = exc_map, 'fig/million-acre-exceedance-map.png', width = 7, height = 3.5)
