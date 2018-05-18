library(tidyverse)
library(rstan)
library(viridis)
library(ggridges)
library(ggrepel)
library(hrbrthemes)
library(patchwork)
library(assertthat)
library(sf)

st_covs <- read_rds('data/processed/st_covs.rds')
cutoff_year <- read_rds('data/processed/cutoff_year.rds')
stan_d <- read_rds('data/processed/stan_d.rds')
min_size <- stan_d$min_size
holdout_burns <- read_rds('data/processed/holdout_burns.rds')
mtbs <- read_rds('data/processed/mtbs.rds')

area_df <- read_rds('data/processed/ecoregions.rds') %>%
  as('Spatial') %>%
  as.data.frame %>%
  tbl_df %>%
  group_by(NA_L3NAME) %>%
  summarize(area = sum(Shape_Area))

# Simulate from the joint predictive distribution -------------------------
# to demonstrate predictions over extremes
# using the ZINB + lognormal model

# we need the predicted counts from the zinb model
# count-preds.rds is generated in count-ppcs.R
zinb_preds <- read_rds('count-preds.rds') %>%
  rename(n_event = value) %>%
  select(-year) %>%
  arrange(iter, id)

# now, having the predicted counts, we need to simulate fire sizes for
# each event, using the predictive distribution from the lognormal model
ln_post <- rstan::extract(read_rds('ba_lognormal_fit.rds'), 
                          pars = c('mu_full', 'scale'))

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


# Inference over total burn area for the test period ----------------------
total_df <- test_preds %>%
  left_join(select(st_covs, id, year)) %>%
  # filter out zero event records (don't contribute to sum)
  filter(n_event > 0, year >= cutoff_year) %>%
  rowwise() %>%
  mutate(total_area = sum(exp(rnorm(n_event, ln_mu, ln_scale)) + min_size)) %>%
  ungroup


actual_totals <- holdout_burns %>%
  select(-geometry) %>%
  as_tibble %>%
  summarize(total_area = sum(R_ACRES),
            total_events = n())

total_df %>%
  group_by(iter) %>%
  summarize(predicted_total_area = sum(total_area),
            predicted_total_events = sum(n_event)) %>%
  ungroup %>%
  mutate(actual_total_area = actual_totals$total_area, 
         actual_total_events = actual_totals$total_events) %>%
  write_csv('data/processed/predicted_totals.csv')
gc()



# Generate derived parameters about the distribution of maxima
nmax_q <- function(p, n, mu, sigma) {
  # quantile function for the n-sample maximum of normally
  # distributed random variables
  erf.inv <- function(x) qnorm((x + 1)/2)/sqrt(2)
  sigma * sqrt(2) * erf.inv(2 * p^(1/n) - 1) + mu
}

# probability of a million acre event, along with medians, 90% pred intervals over block max.
test_preds <- test_preds %>%
  mutate(log_p = n_event * pnorm(log(100000), mean = ln_mu, sd = ln_scale, log.p = TRUE),
         q50 = nmax_q(.5, n = n_event, mu = ln_mu, sigma = ln_scale),
         q95 = nmax_q(.95, n = n_event, mu = ln_mu, sigma = ln_scale),
         q05 = nmax_q(.05, n = n_event, mu = ln_mu, sigma = ln_scale))

holdout_ids <- st_covs$id[st_covs$year >= cutoff_year]
pred_df <- expand.grid(x = 1e6 - min_size,
                       id = holdout_ids,
                       iter = 1:max(test_preds$iter))

iter_max <- 2000
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
n_vs_exc
ggsave(plot = n_vs_exc,
       'fig/number-vs-exceedance.png',
       width = 6, height = 3)

overall_million <- exceedance_df %>%
  # find the cdf for monthly maxima
  mutate(log_p = n_event * pnorm(log(x),
                                 mean = ln_mu,
                                 sd = ln_scale,
                                 log.p = TRUE)) %>%
  left_join(select(st_covs, id, year, rmin)) %>%
  group_by(iter, x) %>%
  # sum to get the joint probability across all time steps
  summarize(total_p = sum(log_p),
            total_events = sum(n_event)) %>%
  ungroup()

p <- ggplot(overall_million, aes(x = 1 - exp(total_p))) +
  geom_histogram(bins = 40, fill = 'darkred', alpha = .7) +
  xlab('Probability of one event > 1,000,000 acres') +
  theme_minimal()
p

# Get the probability of a million+ acre fire
overall_million %>%
  summarize(median(1 - exp(total_p)),
            quantile(1 - exp(total_p), .025),
            quantile(1 - exp(total_p), .975))



# Evaluate interval coverage ----------------------------------------------
# first, get theoretical quantiles for each spatiotemporal unit
test_preds <- test_preds %>%
  mutate(qlo = exp(nmax_q(.005, n = n_event, mu = ln_mu, sigma = ln_scale)),
         qhi = exp(nmax_q(.995, n = n_event, mu = ln_mu, sigma = ln_scale)))

write_rds(test_preds, path = 'test_preds.rds')

interval_df <- test_preds %>%
  group_by(NA_L3NAME, ym) %>%
  summarize(m_qlo = mean(qlo),
            m_qhi = mean(qhi))

max_df <- mtbs %>%
  select(-geometry) %>%
  as.data.frame() %>%
  as_tibble %>%
  group_by(NA_L3NAME, ym) %>%
  summarize(empirical_max = max(R_ACRES - min_size)) %>%
  ungroup %>%
  left_join(distinct(st_covs, NA_L3NAME, NA_L2NAME)) %>%
  group_by(NA_L2NAME) %>%
  mutate(n_events = n())

most_events_ers <- max_df %>%
  distinct(NA_L2NAME, n_events) %>%
  arrange(-n_events)

n_er_to_plot <- 6
er_to_plot <- most_events_ers$NA_L2NAME[1:n_er_to_plot]

interval_df <- interval_df %>%
  left_join(select(st_covs,
                   NA_L3NAME, NA_L2NAME, rmin,
                   ym, year)) %>%
  left_join(max_df) %>%
  filter(year >= cutoff_year) %>%
  mutate(l2_er = tools::toTitleCase(tolower(as.character(NA_L2NAME))),
         l2_er = gsub(' and ', ' & ', l2_er),
         l2_er = gsub('Usa ', '', l2_er),
         l3_er = gsub(' and ', ' & ', NA_L3NAME),
         longer_than_limit = nchar(l3_er) > 32,
         l3_er = substr(l3_er, start = 1, stop = 32),
         l3_er = ifelse(longer_than_limit,
                        paste0(l3_er, '...'),
                        l3_er))

interval_df %>%
  filter(NA_L2NAME %in% er_to_plot) %>%
  ggplot(aes(x = ym, group = NA_L3NAME)) +
  geom_ribbon(aes(ymin = m_qlo, ymax = m_qhi),
              color = NA, alpha = .2,
              fill = 'firebrick') +
  scale_y_log10() +
  theme_minimal() +
  facet_wrap(~ fct_reorder(l2_er, rmin),
             labeller = labeller(.rows = label_wrap_gen(305)),
             nrow = 3) +
  geom_point(aes(y = empirical_max), size = .4) +
  xlab('') +
  ylab('Maximum wildfire size (acres)') +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 7))
ggsave('fig/max-preds-l2-minimal.png', width = 7, height = 4)

# plot for all level 3 ecoregions (for supplement?)
interval_df %>%
  ggplot(aes(x = ym, group = NA_L3NAME)) +
  geom_ribbon(aes(ymin = m_qlo, ymax = m_qhi),
              color = NA,
              fill = 'firebrick', alpha = .6) +
  scale_y_log10() +
  theme_minimal() +
  facet_wrap(~ fct_reorder(l3_er, rmin),
             labeller = labeller(.rows = label_wrap_gen(23))) +
  geom_point(aes(y = empirical_max), size = .3) +
  xlab('') +
  ylab('Maximum event size') +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(size = 5.5),
        axis.text.x = element_text(angle = 90, size = 8),
        axis.text.y = element_text(size = 8))
ggsave('fig/max-preds-l3-all.png', width = 10, height = 10)

# get overall interval coverage stats for block maxima
interval_df %>%
  filter(!is.na(empirical_max)) %>%
  mutate(in_interval = m_qlo <= empirical_max & m_qhi >= empirical_max) %>%
  ungroup %>%
  summarize(p_coverage = mean(in_interval),
            p_max_was_bigger = mean(empirical_max > m_qhi),
            n_max_was_bigger = sum(empirical_max > m_qhi))
