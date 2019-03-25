
# Model comparisons for counts --------------------------------------------
# using the variational approximation
library(cowplot)
library(tidyverse)
library(patchwork)
library(ggrepel)

st_covs <- read_rds('data/processed/st_covs.rds')
cutoff_year <- read_rds('data/processed/cutoff_year.rds')
train_counts <- read_rds('data/processed/train_counts.rds')
stan_d <- read_rds('data/processed/stan_d.rds')

model_fits <- list.files(pattern = '*.fit.*\\.rds')
count_fits <- grep(model_fits, pattern = 'ba_', value = TRUE, invert = TRUE)

# subset to variational fits
count_fits <- grep(count_fits, pattern = 'full', value = TRUE, invert = TRUE)


# data frame for plotting colors
cols <- c('Poisson' = 'green3',
          'Negative binomial' = 'skyblue',
          'ZI Poisson' = 'orange',
          'ZI Negative binomial' = 'hotpink1')

# for each model fit, produce a vector of the holdout log likelihoods
holdout_c_loglik <- list()
train_c_loglik <- list()
train_c_rep <- list()
holdout_c_rep <- list()

for (i in seq_along(count_fits)) {
  post <- rstan::extract(read_rds(count_fits[i]))
  holdout_c_loglik[[i]] <- post$holdout_loglik_c %>%
    reshape2::melt(varnames = c('iter', 'idx')) %>%
    as_tibble %>%
    group_by(iter) %>%
    summarize(value = sum(value)) %>%
    mutate(model = count_fits[i])
  train_c_loglik[[i]] <- post$train_loglik_c %>%
    reshape2::melt(varnames = c('iter', 'idx')) %>%
    as_tibble %>%
    group_by(iter) %>%
    summarize(value = sum(value)) %>%
    mutate(model = count_fits[i])
  pred_counts <- post$count_pred %>%
    reshape2::melt(varnames = c('iter', 'id')) %>%
    as_tibble %>%
    mutate(model = count_fits[i]) %>%
    full_join(dplyr::select(st_covs, id, NA_L3NAME, year, ym))
  train_c_rep[[i]] <- pred_counts %>%
    filter(year < cutoff_year)
  holdout_c_rep[[i]] <- pred_counts %>%
    filter(year >= cutoff_year)
  rm(pred_counts)
  rm(post)
  gc()
}

holdout_c_loglik <- bind_rows(holdout_c_loglik) %>%
  mutate(train = FALSE)

train_c_loglik <- bind_rows(train_c_loglik) %>%
  mutate(train = TRUE)

ll_d <- holdout_c_loglik %>%
  full_join(train_c_loglik) %>%
  mutate(train = ifelse(train == TRUE, 'train', 'test'),
         Distribution = case_when(
           grepl('^nb', .$model) ~ 'Negative binomial',
           grepl('pois', .$model) ~ 'Poisson',
           grepl('zip', .$model) ~ 'ZI Poisson',
           grepl('zinb', .$model) ~ 'ZI Negative binomial')) %>%
  spread(train, value)

loglik_df <- ll_d %>%
  group_by(Distribution) %>%
  summarize(mean_test = median(test),
            sd_test = sd(test)) %>%
  arrange(-mean_test) %>%
  mutate(pretty = paste0(format(mean_test, digits = 0, scientific = FALSE), 
                         ' (', trimws(format(sd_test, digits = 0, scientific = FALSE)), ')')) %>%
  rename(`Holdout log likelihood` = pretty, 
         Model = Distribution) %>%
  select(-ends_with('test')) 

loglik_df %>%
  write_csv('data/processed/count-loglik.csv')


## Generate plots to evaluate distributional assumptions
train_c_rep <- train_c_rep %>%
  bind_rows

holdout_c_rep <- holdout_c_rep %>%
  bind_rows

pr_df <- train_c_rep %>%
  group_by(iter, value, model) %>%
  summarize(n_value = n()) %>%
  ungroup %>%
  group_by(iter, model) %>%
  mutate(total_vals = sum(n_value)) %>%
  ungroup %>%
  group_by(iter, value, model) %>%
  summarize(pr_value = n_value / total_vals)

emp_pr <- train_counts %>%
  mutate(total_vals = n()) %>%
  group_by(n_fire) %>%
  summarize(n_value = n()) %>%
  ungroup %>%
  mutate(pr_value= n_value / sum(n_value))


den_plot <- pr_df %>%
  ungroup %>%
  mutate(Distribution = case_when(
    grepl('^nb', .$model) ~ 'Negative binomial',
    grepl('pois', .$model) ~ 'Poisson',
    grepl('zip', .$model) ~ 'ZI Poisson',
    grepl('zinb', .$model) ~ 'ZI Negative binomial'), 
    Distribution = factor(Distribution, levels = c('ZI Negative binomial', 
                                                   'Negative binomial', 
                                                   'Poisson', 
                                                   'ZI Poisson'))) %>%
  ggplot(aes(value, pr_value, group = iter,
             color = Distribution)) +
  geom_line(alpha = .1) +
  facet_wrap(~model) +
  scale_x_log10(breaks = c(10, 100, 1000)) +
  scale_y_log10(breaks = 10^(c(-5:1))) +
  theme_minimal() +
  scale_color_manual('Count distribution', values = cols) +
  facet_wrap(~Distribution, nrow = 1) +
  geom_point(data = emp_pr,
            aes(x = n_fire, y = pr_value),
            inherit.aes = FALSE, size = .1) +
  theme(legend.position = 'none',
        panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 7)) +
  xlab('Counts > 405 hectares') +
  ylab('Probability')

# proportion of zero observations
ppc_counts <- train_c_rep %>%
  group_by(iter, model) %>%
  summarize(p_zero = mean(value == 0),
            max_count = max(value),
            sum_count = sum(value))

test_ppc <- holdout_c_rep %>%
  group_by(iter, model) %>%
  summarize(p_zero = mean(value == 0),
            max_count = max(value),
            sum_count = sum(value))

zero_plot <- ppc_counts %>%
  mutate(train = 'train') %>%
  ungroup %>%
  full_join(test_ppc %>% mutate(train = 'test')) %>%
  mutate(Distribution = case_when(
    grepl('^nb', .$model) ~ 'Negative binomial',
    grepl('pois', .$model) ~ 'Poisson',
    grepl('zip', .$model) ~ 'ZI Poisson',
    grepl('zinb', .$model) ~ 'ZI Negative binomial'), 
    Distribution = factor(Distribution, levels = c('ZI Negative binomial', 
                                                   'Negative binomial', 
                                                   'Poisson', 
                                                   'ZI Poisson'))) %>%
  dplyr::select(-max_count, -sum_count, -model) %>%
  spread(train, p_zero) %>%
  ggplot(aes(x = train, y = test, color = Distribution)) +
  facet_wrap(~ Distribution, nrow = 1) +
  theme_minimal() +
  geom_point(alpha = .4) +
  scale_fill_manual('Count distribution', values = cols) +
  geom_vline(xintercept = mean(stan_d$counts == 0),
             linetype = 'dashed', color = 'black') +
  geom_hline(yintercept = mean(stan_d$holdout_c == 0),
             linetype = 'dashed', color = 'black') +
  xlab('Pr(0): train') +
  ylab('Pr(0): test') +
  scale_color_manual('Count distribution', values = cols) +
  theme(legend.position = 'none') +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.text = element_text(size = 7)
  )

max_plot <- ppc_counts %>%
  mutate(train = 'train') %>%
  ungroup %>%
  full_join(test_ppc %>% mutate(train = 'test')) %>%
  mutate(Distribution = case_when(
    grepl('^nb', .$model) ~ 'Negative binomial',
    grepl('pois', .$model) ~ 'Poisson',
    grepl('zip', .$model) ~ 'ZI Poisson',
    grepl('zinb', .$model) ~ 'ZI Negative binomial'), 
    Distribution = factor(Distribution, levels = c('ZI Negative binomial', 
                                                   'Negative binomial', 
                                                   'Poisson', 
                                                   'ZI Poisson'))) %>%
  dplyr::select(-p_zero, -sum_count, -model) %>%
  spread(train, max_count) %>%
  ggplot(aes(x = train, y = test, color = Distribution)) +
  theme_minimal() +
  geom_point(alpha = .4) +
  facet_wrap(~ Distribution, nrow = 1) +
  scale_fill_manual('Count distribution', values = cols) +
  xlab('max: train') +
  ylab('max: test') +
  scale_color_manual('Count distribution', values = cols) +
  # scale_x_log10() +
  # scale_y_log10() +
  geom_vline(xintercept = max(stan_d$counts),
             linetype = 'dashed', color = 'black') +
  geom_hline(yintercept = max(stan_d$holdout_c),
             linetype = 'dashed', color = 'black') +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 7))

sum_plot <- ppc_counts %>%
  mutate(train = 'train') %>%
  ungroup %>%
  full_join(test_ppc %>% mutate(train = 'test')) %>%
  mutate(Distribution = case_when(
    grepl('^nb', .$model) ~ 'Negative binomial',
    grepl('pois', .$model) ~ 'Poisson',
    grepl('zip', .$model) ~ 'ZI Poisson',
    grepl('zinb', .$model) ~ 'ZI Negative binomial'), 
    Distribution = factor(Distribution, levels = c('ZI Negative binomial', 
                                                   'Negative binomial', 
                                                   'Poisson', 
                                                   'ZI Poisson'))) %>%
  dplyr::select(-max_count, -p_zero, -model) %>%
  spread(train, sum_count) %>%
  ggplot(aes(x = train, y = test, color = Distribution)) +
  theme_minimal() +
  facet_wrap(~ Distribution, nrow = 1) +
  geom_point(alpha = .4) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  scale_fill_manual('Count distribution', values = cols) +
  xlab('sum: train') +
  ylab('sum: test') +
  scale_color_manual('Count distribution', values = cols) +
  geom_vline(xintercept = sum(stan_d$counts),
             linetype = 'dashed', color = 'black') +
  geom_hline(yintercept = sum(stan_d$holdout_c),
             linetype = 'dashed', color = 'black') +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 7))

cowplot::plot_grid(den_plot, zero_plot, max_plot, sum_plot, 
                   ncol = 1, rel_heights = c(1.3, 1, 1, 1))
ggsave('fig/figure_4.pdf', width = 6, height = 4.5)
