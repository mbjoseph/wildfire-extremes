library(tidyverse)
library(ggrepel)

st_covs <- read_rds('data/processed/st_covs.rds')
cutoff_year <- read_rds('data/processed/cutoff_year.rds')
train_burns <- read_rds('data/processed/train_burns.rds')
holdout_burns <- read_rds('data/processed/holdout_burns.rds')
stan_d <- read_rds('data/processed/stan_d.rds')

model_fits <- list.files(pattern = '*.fit.*\\.rds')

burn_area_fits <- grep(model_fits, pattern = 'ba_', value = TRUE)

hectares_per_acre <- 0.404686

# for each model fit, produce a vector of the holdout log likelihoods
holdout_ba_loglik <- list()
train_ba_loglik <- list()
train_ba_rep <- list()
holdout_ba_rep <- list()

for (i in seq_along(burn_area_fits)) {
  post <- rstan::extract(read_rds(burn_area_fits[i]))
  holdout_ba_loglik[[i]] <- post$holdout_loglik_b %>%
    reshape2::melt(varnames = c('iter', 'idx')) %>%
    as_tibble %>%
    group_by(iter) %>%
    summarize(value = sum(value)) %>%
    mutate(model = burn_area_fits[i])

  train_ba_loglik[[i]] <- post$loglik_f %>%
    reshape2::melt(varnames = c('iter', 'idx')) %>%
    as_tibble %>%
    group_by(iter) %>%
    summarize(value = sum(value)) %>%
    mutate(model = burn_area_fits[i])

  holdout_ba_rep[[i]] <- post$holdout_rep %>%
    reshape2::melt(varnames = c('iter', 'idx')) %>%
    as_tibble %>%
    mutate(model = burn_area_fits[i])

  train_ba_rep[[i]] <- post$size_rep %>%
    reshape2::melt(varnames = c('iter', 'idx')) %>%
    as_tibble %>%
    mutate(model = burn_area_fits[i])

  rm(post)
  gc()
}

holdout_ba_loglik <- bind_rows(holdout_ba_loglik) %>%
  mutate(train = FALSE)

train_ba_loglik <- bind_rows(train_ba_loglik) %>%
  mutate(train = TRUE)

# data frame for plotting colors
cols <- c('Gamma' = '#66c2a5',
          'Weibull' = '#fc8d62',
          'Generalized Pareto' = '#8da0cb',
          'Tapered Pareto' = '#a6d854',
          'Lognormal' = '#e78ac3')

loglik_d <- holdout_ba_loglik %>%
  full_join(train_ba_loglik) %>%
  mutate(train = ifelse(train == TRUE, 'train', 'test'),
         Distribution = case_when(
           grepl('gamma', .$model) ~ 'Gamma',
           grepl('lognormal', .$model) ~ 'Lognormal',
           grepl('_pareto', .$model) ~ 'Generalized Pareto',
           grepl('_tpareto', .$model) ~ 'Tapered Pareto',
           grepl('weibull', .$model) ~ 'Weibull'
         )) %>%
  spread(train, value)

label_d <- loglik_d %>%
  group_by(Distribution) %>%
  summarize(med_test = median(test),
            med_train = median(train))

loglik_d %>%
  group_by(Distribution) %>%
  summarize(mean_test = median(test),
            sd_test = sd(test)) %>%
  arrange(-mean_test) %>%
  mutate(pretty = paste0(format(mean_test, digits = 0, scientific = FALSE), 
                         ' (', trimws(format(sd_test, digits = 0, scientific = FALSE)), ')')) %>%
  rename(`Holdout log likelihood` = pretty, 
         Model = Distribution) %>%
  select(-ends_with('test')) %>%
  write_csv('data/processed/burn-area-loglik.csv')

rm(holdout_ba_loglik)
rm(train_ba_loglik)
gc()

## Generate plots to evaluate distributional assumptions
train_ba_rep <- train_ba_rep %>%
  bind_rows


my_theme <- theme_minimal() + 
  theme(legend.position = 'none', 
        panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 6), 
        axis.title = element_text(size = 9))

# Whole distribution
den_plot <- train_ba_rep %>%
  filter(iter < 300) %>%
  mutate(Distribution = case_when(
    grepl('gamma', .$model) ~ 'Gamma',
    grepl('lognormal', .$model) ~ 'Lognormal',
    grepl('_pareto', .$model) ~ 'Generalized Pareto',
    grepl('_tpareto', .$model) ~ 'Tapered Pareto',
    grepl('weibull', .$model) ~ 'Weibull'), 
    Distribution = factor(Distribution, levels = c('Lognormal', 
                                                   'Generalized Pareto', 
                                                   'Tapered Pareto', 
                                                   'Weibull', 
                                                   'Gamma'))) %>%
  ggplot(aes(x = value * hectares_per_acre,
             color = Distribution,
             group = interaction(Distribution, iter))) +
  my_theme +
  scale_x_log10() +
  stat_density(aes(color = Distribution), geom="line", alpha = .1, position = 'identity') +
  scale_color_manual('Burn area distribution', values = cols) +
  facet_wrap(~Distribution, nrow = 1) +
  stat_density(geom = 'line', position = 'identity',
               inherit.aes = FALSE,
               aes(x = actual_sizes * hectares_per_acre),
               data = tibble(actual_sizes = stan_d$sizes)) +
  coord_cartesian(xlim = c(min(stan_d$sizes), 700000)) +
  geom_rug(inherit.aes = FALSE,
           aes(x = actual_sizes),
           data = tibble(actual_sizes = stan_d$sizes),
           alpha = .1) +
  xlab('Burn area exceedance (hectares)') +
  ylab('Density')

# "Survival" plots for top 3 models
obs_ccdf <- tibble(values = sort(stan_d$sizes)) %>%
  mutate(ccdf = seq(n(), 1, -1) / n()) %>%
  filter(ccdf <= .01)
gc()

tail_plot <- train_ba_rep %>%
  filter(iter < 300) %>%
  group_by(iter, model) %>%
  arrange(iter, model, value) %>%
  mutate(ccdf = seq(n(), 1, -1) / n()) %>%
  filter(ccdf <= .01) %>%
  mutate(value = value * hectares_per_acre) %>%
  group_by(ccdf, model) %>%
  summarize(lo = quantile(value, .025),
            q25 = quantile(value, .25),
            med = median(value),
            q75 = quantile(value, .75),
            hi = quantile(value, .975)) %>%
  ungroup %>%
  mutate(Distribution = case_when(
    grepl('gamma', .$model) ~ 'Gamma',
    grepl('lognormal', .$model) ~ 'Lognormal',
    grepl('_pareto', .$model) ~ 'Generalized Pareto',
    grepl('_tpareto', .$model) ~ 'Tapered Pareto',
    grepl('weibull', .$model) ~ 'Weibull'
  ), 
  Distribution = factor(Distribution, levels = c('Lognormal', 
                                                 'Generalized Pareto', 
                                                 'Tapered Pareto', 
                                                 'Weibull', 
                                                 'Gamma'))) %>%
  ggplot(aes(med, ccdf, color = Distribution)) +
  my_theme +
  scale_color_manual('Burn area distribution', values = cols) +
  geom_point(shape = 1) +
  geom_errorbarh(aes(xmin = lo, xmax = hi)) +
  geom_errorbarh(aes(xmin = q25, xmax = q75), size = 1.25) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~Distribution, nrow = 1) +
  coord_cartesian(xlim = c(1e5 * hectares_per_acre, 1e8 * hectares_per_acre),
                  ylim = c(1e-4, .01)) +
  geom_point(data = obs_ccdf, color = 'black', aes(x = values * hectares_per_acre), shape = 1, size = .5) +
  xlab('Burn area exceedance (hectares)') +
  ylab('CCDF') +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )


# Maxima and cumulative burn area -----------------------------------------
train_max <- train_ba_rep %>%
  group_by(iter, model) %>%
  summarize(max_pred = max(value),
            cum_pred = sum(value)) %>%
  mutate(train = TRUE)

holdout_max <- holdout_ba_rep %>%
  bind_rows  %>%
  group_by(iter, model) %>%
  summarize(max_pred = max(value),
            cum_pred = sum(value)) %>%
  mutate(train = FALSE)

gc()

actual_maxima <- tibble(train = max(stan_d$sizes),
                        test = max(stan_d$holdout_b))

actual_sums <- tibble(train = sum(stan_d$sizes),
                      test = sum(stan_d$holdout_b))

max_plot <- train_max %>%
  full_join(holdout_max) %>%
  dplyr::select(-cum_pred) %>%
  mutate(train = ifelse(train, 'train', 'test')) %>%
  spread(train, max_pred) %>%
  ungroup() %>%
  mutate(Distribution = case_when(
    grepl('gamma', .$model) ~ 'Gamma',
    grepl('lognormal', .$model) ~ 'Lognormal',
    grepl('_pareto', .$model) ~ 'Generalized Pareto',
    grepl('_tpareto', .$model) ~ 'Tapered Pareto',
    grepl('weibull', .$model) ~ 'Weibull'
  ), 
  Distribution = factor(Distribution, levels = c('Lognormal', 
                                                 'Generalized Pareto', 
                                                 'Tapered Pareto', 
                                                 'Weibull', 
                                                 'Gamma'))) %>%
  ggplot(aes(train * hectares_per_acre, test * hectares_per_acre, color = Distribution)) +
  my_theme +
  scale_color_manual('Burn area distribution', values = cols) +
  geom_point(alpha = .4) +
  xlab('max(area): train') +
  ylab('max(area): test') +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~ Distribution, nrow = 1) +
  geom_vline(aes(xintercept = train * hectares_per_acre),
             data = actual_maxima, color = 'black', linetype = 'dashed') +
  geom_hline(aes(yintercept = test * hectares_per_acre),
             data = actual_maxima, color = 'black', linetype = 'dashed') +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

sum_plot <- train_max %>%
  full_join(holdout_max) %>%
  dplyr::select(-max_pred) %>%
  mutate(train = ifelse(train, 'train', 'test')) %>%
  spread(train, cum_pred) %>%
  ungroup() %>%
  mutate(Distribution = case_when(
    grepl('gamma', .$model) ~ 'Gamma',
    grepl('lognormal', .$model) ~ 'Lognormal',
    grepl('_pareto', .$model) ~ 'Generalized Pareto',
    grepl('_tpareto', .$model) ~ 'Tapered Pareto',
    grepl('weibull', .$model) ~ 'Weibull'
  ), 
  Distribution = factor(Distribution, levels = c('Lognormal', 
                                                 'Generalized Pareto', 
                                                 'Tapered Pareto', 
                                                 'Weibull', 
                                                 'Gamma'))) %>%
  ggplot(aes(train * hectares_per_acre, test * hectares_per_acre, color = Distribution)) +
  my_theme +
  scale_color_manual('Burn area distribution', values = cols) +
  geom_point(alpha = .4) +
  xlab('sum(area): train') +
  ylab('sum(area): test') +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~ Distribution, nrow = 1) +
  geom_vline(aes(xintercept = train * hectares_per_acre),
             data = actual_sums, color = 'black', linetype = 'dashed') +
  geom_hline(aes(yintercept = test * hectares_per_acre),
             data = actual_sums, color = 'black', linetype = 'dashed') +
  theme(strip.background = element_blank(),
    strip.text.x = element_blank())

cowplot::plot_grid(den_plot, tail_plot, max_plot, sum_plot, 
                   nrow = 4, rel_heights = c(1.3, 1, 1, 1))
ggsave(filename = 'fig/figure_7.pdf', width = 6, height = 4.25)





# Evaluate posterior predictive interval coverage in test set -------------
test_pred_intervals <- holdout_ba_rep %>%
  bind_rows %>%
  group_by(idx, model) %>%
  summarize(med = median(value),
            lo = quantile(value, .025),
            hi = quantile(value, .975)) %>%
  mutate(true_exceedance = stan_d$holdout_b[idx],
         true_in_interval = lo <= true_exceedance & hi >= true_exceedance,
         NA_L3NAME = holdout_burns$NA_L3NAME[idx])

write_csv(test_pred_intervals, 'data/processed/area_coverage.csv')
