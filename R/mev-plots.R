library(lubridate)
library(tidyverse)
library(patchwork)
library(assertthat)
library(sf)

st_covs <- read_rds('data/processed/st_covs.rds')
cutoff_year <- read_rds('data/processed/cutoff_year.rds')
stan_d <- read_rds('data/processed/stan_d.rds')
min_size <- stan_d$min_size
holdout_burns <- read_rds('data/processed/holdout_burns.rds')
mtbs <- read_rds('data/processed/mtbs.rds')


# Simulate from the joint predictive distribution -------------------------
# to demonstrate predictions over extremes
# using the ZINB + lognormal model
zinb_preds <- read_rds('count-preds.rds') %>%
  rename(n_event = value) %>%
  select(-year) %>%
  arrange(iter, id)

ln_post <- rstan::extract(read_rds('ba_lognormal_fit.rds'), 
                          pars = c('mu_full', 'scale'))

ln_mu <- ln_post$mu_full %>%
  reshape2::melt(varnames = c('iter', 'id'), value.name = 'ln_mu') %>%
  arrange(iter, id)

ln_sigma <- ln_post$scale %>%
  reshape2::melt(varnames = c('iter'), value.name = 'ln_scale')

# merge lognormal mean and scale draws with counts
# (bind cols is much faster than *_join)
test_preds <- bind_cols(zinb_preds, ln_mu) %>%
  left_join(ln_sigma)

# ensure that indices match for merge
assert_that(all(test_preds$iter == test_preds$iter1))
assert_that(all(test_preds$id == test_preds$id1))

# remove unneeded columns
test_preds <- test_preds %>%
  select(-ends_with('1'))

# subset to test period
test_preds <- test_preds %>%
  left_join(select(st_covs, id, year)) %>%
  filter(year >= cutoff_year)
gc()

# Inference over total burn area for the test period ----------------------
total_df <- test_preds %>%
  # filter out zero event records (don't contribute to sum)
  filter(n_event > 0) %>%
  rowwise() %>%
  mutate(total_area = sum(exp(rnorm(n_event, ln_mu, ln_scale)) + min_size)) %>%
  ungroup

actual_totals <- holdout_burns %>%
  select(-geometry) %>%
  as_tibble %>%
  summarize(total_area = sum(Acres),
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


# Simulating maxima from posterior predictive dist ------------------------
pr_zero <- test_preds %>%
  group_by(ym, NA_L3NAME) %>%
  summarize(pr_zero = mean(n_event == 0))


max_df <- test_preds %>%
  # filter out zero event records (don't contribute to max)
  filter(n_event > 0) %>%
  rowwise() %>%
  mutate(max_size = max(exp(rnorm(n_event, ln_mu, ln_scale)) + min_size), 
         pr_over_million = mean(exp(rnorm(n_event, ln_mu, ln_scale)) + min_size >= 1e6)) %>%
  ungroup

zero_df <- test_preds %>%
  filter(n_event == 0) %>%
  mutate(pr_over_million = 0, 
         max_size = min_size)

mill_df <- max_df %>%
  full_join(zero_df)

mill_summ <- mill_df %>%
  group_by(NA_L3NAME, iter) %>% # average over timesteps
  summarize(any_million_acre_events = any(max_size > 1e6)) %>%
  ungroup %>%
  group_by(NA_L3NAME) %>%
  summarize(pr_million_acre_event = mean(any_million_acre_events))
write_rds(mill_summ, path = 'data/processed/mill_summ.rds')

overall_summ <- mill_df %>%
  group_by(iter) %>% # any million acre events in any location, any timestep
  summarize(any_million_acre_events = any(max_size > 1e6)) %>%
  ungroup %>%
  summarize(pr_million_acre_event = mean(any_million_acre_events))


# Generate derived parameters about the distribution of maxima ----------------
nmax_q <- function(p, n, mu, sigma) {
  # quantile function for the n-sample maximum of normally
  # distributed random variables
  erf.inv <- function(x) qnorm((x + 1)/2)/sqrt(2)
  sigma * sqrt(2) * erf.inv(2 * p^(1/n) - 1) + mu
}


# just a quick check to make sure this is implemented correctly:
# points should fall on or very near one to one line
n_vec <- c(1, 5, 10, 20, 100, 1000)
ps <- c(.5, .9)
expand.grid(n = n_vec, iter = 1:10000) %>%
  as_tibble %>%
  rowwise %>%
  mutate(ymax = max(rnorm(n))) %>%
  ungroup %>%
  group_by(n) %>%
  summarize(q5 = median(ymax), 
            q1 = quantile(ymax, .1), 
            q9 = quantile(ymax, .9)) %>%
  gather(quantile, value, -n) %>%
  mutate(quantile = case_when(
    .$quantile == "q5" ~ .5, 
    .$quantile == "q1" ~ .1, 
    .$quantile == "q9" ~ .9
  ), 
  analytic_value = nmax_q(quantile, n, 0, 1)) %>%
  ggplot(aes(x = analytic_value, y = value)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed') + 
  xlab("Analytic n-sample maximum quantile") + 
  ylab("Empirical n-sample maximum quantile")


# log_p: CDF of maximum distribution, evaluated at an exceedance of one million - 1e3
# (the probability that the maximum is less than one million acres)
test_preds <- test_preds %>%
  mutate(log_p = n_event * pnorm(log(1e6 - 1e3), mean = ln_mu, sd = ln_scale, log.p = TRUE))

holdout_ids <- st_covs$id[st_covs$year >= cutoff_year]
pred_df <- expand.grid(x = 1e6 - min_size,
                       id = holdout_ids,
                       iter = 1:max(test_preds$iter))

exceedance_df <- test_preds %>%
  filter(id %in% holdout_ids) %>%
  full_join(pred_df)
 
overall_million <- exceedance_df %>%
  left_join(select(st_covs, id, year, rmin)) %>%
  group_by(iter, x) %>%
  # sum over ecoregions and timesteps to get the joint probability that all 
  # maxima are less than one million across all time steps
  summarize(total_p = sum(log_p)) %>%
  ungroup()

# probability of a million+ acre fire is the complement of the log probability
# that all block maxima are less than one million
overall_million %>%
  summarize(med = median(1 - exp(total_p)),
            lo = quantile(1 - exp(total_p), .025),
            hi = quantile(1 - exp(total_p), .975)) %>%
  write_csv('data/processed/million-acre-prob.csv')

# probabilities of million acre fires for each ecoregion month
million_er_mon <- exceedance_df %>%
  group_by(x, NA_L3NAME, ym) %>%
  # summarize posterior probability of million acre fires over post. draws
  summarize(med = median(1 - exp(log_p)),
            lo = quantile(1 - exp(log_p), .025),
            hi = quantile(1 - exp(log_p), .975), 
            expected_n = mean(n_event)) %>%
  ungroup() %>%
  left_join(select(st_covs, NA_L3NAME, NA_L2NAME, NA_L1NAME, ym, rmin, month, year))

million_er_mon %>%
  write_csv('data/processed/million-er-mon.csv')

# Million acre probability by ym
ym_million_prob <- bind_cols(zinb_preds, ln_mu) %>%
  left_join(ln_sigma) %>%
  select(-ends_with('1')) %>%
  left_join(select(st_covs, id, year)) %>%
  mutate(log_p = n_event * pnorm(log(1e6 - 1e3), mean = ln_mu, sd = ln_scale, log.p = TRUE)) %>%
  group_by(iter, ym) %>%
  summarize(total_p = sum(log_p)) %>%
  mutate(pr_million = 1 - exp(total_p))

ym_million_prob %>%
  group_by(ym) %>%
  summarize(med = median(pr_million), 
            lo = quantile(pr_million, .1), 
            hi = quantile(pr_million, .9)) %>%
  mutate(test_set = ifelse(ym >= 2010, 'test', 'train')) %>%
  ggplot(aes(x = ym, y = med, fill = test_set)) + 
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .5) +
  geom_line(rel_min_height = .01) + 
  theme_minimal() + 
  ylab("Million acre wildfire probability") + 
  xlab("") + 
  scale_fill_manual(values = c('red', 'black')) + 
  theme(panel.grid.minor = element_blank(), 
        legend.position = 'none') + 
  annotate(geom = 'text', x = 2007, y = .15, label = 'Train') + 
  annotate(geom = 'text', x = 2013, y = .15, label = 'Test', color = 'red')
ggsave('fig/million-acre-probs.png', width = 4, height = 3)


# Evaluate interval coverage ----------------------------------------------
# first, get theoretical quantiles for each spatiotemporal unit
test_preds <- test_preds %>%
  mutate(qlo = min_size + exp(nmax_q(.025, n = n_event, mu = ln_mu, sigma = ln_scale)),
         qhi = min_size + exp(nmax_q(.975, n = n_event, mu = ln_mu, sigma = ln_scale)), 
         qvlo = min_size + exp(nmax_q(.005, n = n_event, mu = ln_mu, sigma = ln_scale)),
         qvhi = min_size + exp(nmax_q(.995, n = n_event, mu = ln_mu, sigma = ln_scale)), 
         q25 = min_size + exp(nmax_q(.25, n = n_event, mu = ln_mu, sigma = ln_scale)), 
         q75 = min_size + exp(nmax_q(.75, n = n_event, mu = ln_mu, sigma = ln_scale)), 
         q50 = min_size + exp(nmax_q(.5, n = n_event, mu = ln_mu, sigma = ln_scale)))
write_rds(test_preds, path = 'test_preds.rds')




# Visualize extremes predictions for some ecoregions ----------------------
interval_df <- test_preds %>%
  group_by(NA_L3NAME, ym) %>%
  summarize(m_qlo = mean(qlo),
           m_qhi = mean(qhi), 
           m_med = mean(q50),
           m_qvlo = mean(qvlo), 
           m_qvhi = mean(qvhi))

max_df <- mtbs %>%
  select(-geometry) %>%
  as.data.frame() %>%
  as_tibble %>%
  group_by(NA_L3NAME, ym) %>%
  mutate(empirical_max = max(Acres), 
         is_empirical_max = Acres == empirical_max) %>%
  ungroup %>%
  filter(is_empirical_max) %>%
  left_join(distinct(st_covs, NA_L3NAME, NA_L2NAME)) %>%
  group_by(NA_L2NAME) %>%
  mutate(n_events = n())

interval_df <- interval_df %>%
  left_join(select(st_covs,
                   NA_L3NAME, NA_L2NAME, rmin, month,
                   ym, year)) %>%
  left_join(max_df) %>%
  filter(year >= cutoff_year) %>%
  mutate(l2_er = str_to_title(NA_L2NAME),
         l2_er = gsub(' and ', ' & ', l2_er),
         l2_er = gsub('Usa ', '', l2_er),
         l3_er = gsub(' and ', ' & ', NA_L3NAME),
         longer_than_limit = nchar(l3_er) > 32,
         l3_er = substr(l3_er, start = 1, stop = 32),
         l3_er = ifelse(longer_than_limit,
                        paste0(l3_er, '...'),
                        l3_er))

# plot for level 3 ecoregions
hectares_per_acre <- 0.404686
interval_df %>%
  group_by(NA_L3NAME) %>%
  mutate(total_n = sum(!is.na(empirical_max))) %>%
  filter(total_n > 26) %>%
  ggplot(aes(x = ym, group = NA_L3NAME)) +
  geom_ribbon(aes(ymin = m_qlo * hectares_per_acre, ymax = m_qhi * hectares_per_acre),
              color = NA,
              fill = 'firebrick', alpha = .7) +
  geom_ribbon(aes(ymin = m_qvlo * hectares_per_acre, ymax = m_qvhi * hectares_per_acre),
              color = NA,
              fill = 'firebrick', alpha = .3) +
  scale_y_log10() +
  theme_minimal() +
  facet_wrap(~ fct_reorder(l3_er, rmin),
             labeller = labeller(.rows = label_wrap_gen(23))) +
  geom_point(aes(y = empirical_max * hectares_per_acre), size = .5) +
  xlab('') +
  ylab('Maximum fire size (hectares)') +
  theme(panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 90))
ggsave('fig/max-preds-l3.png', width = 7, height = 4)

# get overall interval coverage stats for block maxima
write_csv(interval_df, 'data/processed/mev_intervals.csv')
