
# Model comparisons for counts --------------------------------------------
library(tidyverse)
library(cowplot)
library(extraDistr)
library(ggridges)
library(ggExtra)
library(patchwork)

source('R/02-explore.R')
source('R/make-stan-d.R')


# load posterior predictions for counts
train_c_rep <- read_rds('train_c_rep.rds')
holdout_c_rep <- read_rds('holdout_c_rep.rds')

# subset for memory management
train_c_rep <- filter(train_c_rep, iter < 500)
holdout_c_rep <- filter(holdout_c_rep, iter < 500)

# Total burn areas by year
train_y_rep <- train_c_rep %>%
  group_by(year, iter, model) %>%
  summarize(total_pred = sum(value),
            max_pred = max(value),
            p_zero = mean(value == 0)) %>%
  ungroup %>%
  mutate(train = 'train')

holdout_y_rep <- holdout_c_rep %>%
  group_by(year, iter, model) %>%
  summarize(total_pred = sum(value),
            max_pred = max(value),
            p_zero = mean(value == 0)) %>%
  ungroup %>%
  mutate(train = 'test')

actual_vals <- count_df %>%
  group_by(year) %>%
  summarize(max_pred = max(n_fire),
            p_zero = mean(n_fire == 0),
            total_pred = sum(n_fire)) %>%
  gather(target, value, -year)

train_y_rep %>%
  full_join(holdout_y_rep) %>%
  filter(model == 'nb_fit.rds') %>%
  gather(target, value, -year, -iter, -train, -model) %>%
  ggplot(aes(x = year, y = value, color = train)) +
  geom_jitter(alpha = .1, width = .4) +
  facet_wrap(~target, scales = 'free') +
  scale_y_log10() +
  geom_point(data = actual_vals, color = 'black')
ggsave('fig/count-ppcs-by-year.png')



## Inspect by ecoregion and year for holdout data
train_full_rep <- train_c_rep %>%
  group_by(ym, model, NA_L3NAME) %>%
  summarize(med = median(value),
            lo = quantile(value, .025),
            hi = quantile(value, .975)) %>%
  ungroup

holdout_full_rep <- holdout_c_rep %>%
  group_by(ym, model, NA_L3NAME) %>%
  summarize(med = median(value),
            lo = quantile(value, .025),
            hi = quantile(value, .975)) %>%
  ungroup

holdout_full_rep %>%
  ggplot(aes(ym, med)) +
  geom_line() +
  facet_wrap(~NA_L3NAME) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .4) +
  scale_y_log10(breaks = c(0, 1, 10, 100, 1000, 10000)) +
  theme_minimal() +
  geom_point(data = holdout_counts, aes(y = n_fire), col = 'red', size = .5)
ggsave('fig/count-preds.png', width = 18, height = 10)

# zoom in on problematic ecoregions
train_full_rep %>%
  mutate(train = 'train') %>%
  full_join(holdout_full_rep %>% mutate(train = 'test')) %>%
  full_join(er_df) %>%
  ggplot(aes(ym, med)) +
  geom_line() +
  facet_wrap(~NA_L3NAME) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = train), alpha = .8) +
  scale_fill_manual(values = c('black', 'darkblue')) +
  scale_y_log10(breaks = c(0, 1, 10, 100, 1000, 10000)) +
  theme_minimal() +
  geom_point(data = count_df, aes(y = n_fire), col = 'red', size = .5)
ggsave('fig/count-preds-full.png', width = 30, height = 10)

