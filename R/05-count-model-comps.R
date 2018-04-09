
# Model comparisons for counts --------------------------------------------
# using the variational approximation
library(tidyverse)
library(cowplot)
library(extraDistr)
library(ggridges)
library(ggExtra)
library(patchwork)
library(ggrepel)

source('R/02-explore.R')
source('R/make-stan-d.R')

# possibly fetch from Amazon
#system("aws s3 cp s3://earthlab-mjoseph .. --recursive --exclude '*' --include '*.rds'")

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

ll_labels <- ll_d %>%
  group_by(Distribution) %>%
  summarize(test = median(test),
            train = median(train))

ll_d %>%
  ggplot(aes(x = train, y = test, color = Distribution)) +
  theme_minimal() +
  geom_point(alpha = .5) +
  xlab('Log likelihood: training set') +
  ylab('Log likelihood: test set') +
  scale_color_manual('Count distribution', values = cols) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(panel.grid.minor = element_blank(),
        legend.position = 'none') +
  geom_text_repel(data = ll_labels,
                  aes(label = Distribution),
                  color = 'black')
ggsave(filename = 'fig/loglik-counts.png', width = 6, height = 4)



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
    grepl('zinb', .$model) ~ 'ZI Negative binomial')) %>%
  ggplot(aes(value, pr_value, group = iter,
             color = Distribution)) +
  geom_line(alpha = .1) +
  facet_wrap(~model) +
  scale_x_log10() +
  scale_y_log10(breaks = 10^(c(-5:1))) +
  theme_minimal() +
  scale_color_manual('Count distribution', values = cols) +
  facet_wrap(~Distribution, nrow = 1) +
  geom_point(data = emp_pr,
            aes(x = n_fire, y = pr_value),
            inherit.aes = FALSE) +
  theme(legend.position = 'none',
        panel.grid.minor = element_blank()) +
  xlab('Counts > 1000 acres') +
  ylab('Pr(counts > 1000 acres)')

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
    grepl('zinb', .$model) ~ 'ZI Negative binomial')) %>%
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
  xlab('Pr(0): training data') +
  ylab('Pr(0): test data') +
  scale_color_manual('Count distribution', values = cols) +
  theme(legend.position = 'none') +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.grid.minor = element_blank()
  )

max_plot <- ppc_counts %>%
  mutate(train = 'train') %>%
  ungroup %>%
  full_join(test_ppc %>% mutate(train = 'test')) %>%
  mutate(Distribution = case_when(
    grepl('^nb', .$model) ~ 'Negative binomial',
    grepl('pois', .$model) ~ 'Poisson',
    grepl('zip', .$model) ~ 'ZI Poisson',
    grepl('zinb', .$model) ~ 'ZI Negative binomial')) %>%
  dplyr::select(-p_zero, -sum_count, -model) %>%
  spread(train, max_count) %>%
  ggplot(aes(x = train, y = test, color = Distribution)) +
  theme_minimal() +
  geom_point(alpha = .4) +
  facet_wrap(~ Distribution, nrow = 1) +
  scale_fill_manual('Count distribution', values = cols) +
  xlab('max(count): training data') +
  ylab('max(count): test data') +
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
        panel.grid.minor = element_blank())

sum_plot <- ppc_counts %>%
  mutate(train = 'train') %>%
  ungroup %>%
  full_join(test_ppc %>% mutate(train = 'test')) %>%
  mutate(Distribution = case_when(
    grepl('^nb', .$model) ~ 'Negative binomial',
    grepl('pois', .$model) ~ 'Poisson',
    grepl('zip', .$model) ~ 'ZI Poisson',
    grepl('zinb', .$model) ~ 'ZI Negative binomial')) %>%
  dplyr::select(-max_count, -p_zero, -model) %>%
  spread(train, sum_count) %>%
  ggplot(aes(x = train, y = test, color = Distribution)) +
  theme_minimal() +
  facet_wrap(~ Distribution, nrow = 1) +
  geom_point(alpha = .4) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  scale_fill_manual('Count distribution', values = cols) +
  xlab('sum(count): training data') +
  ylab('sum(count): test data') +
  scale_color_manual('Count distribution', values = cols) +
  geom_vline(xintercept = sum(stan_d$counts),
             linetype = 'dashed', color = 'black') +
  geom_hline(yintercept = sum(stan_d$holdout_c),
             linetype = 'dashed', color = 'black') +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.minor = element_blank())

plot_grid(den_plot, zero_plot, max_plot, sum_plot, ncol = 1)
ggsave('fig/ppc-counts.png', width = 9, height = 7)
