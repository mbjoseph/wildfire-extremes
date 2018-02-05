
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

# Total counts by year
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
  filter(model == 'zinb_fit.rds') %>%
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
  filter(model == 'zinb_fit.rds') %>%
  left_join(holdout_counts) %>%
  mutate(in_interval = n_fire >= lo & n_fire <= hi,
         l3_er = gsub(' and ', ' & ', NA_L3NAME),
         longer_than_limit = nchar(l3_er) > 32,
         l3_er = substr(l3_er, start = 1, stop = 32),
         l3_er = ifelse(longer_than_limit,
                        paste0(l3_er, '...'),
                        l3_er)) %>%
  ggplot(aes(ym, med, color = in_interval)) +
  facet_wrap(~reorder(l3_er, in_interval),
             labeller = labeller(.rows = label_wrap_gen(25)),
             nrow = 12) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .2, color = NA) +
  scale_y_log10(breaks = c(0, 1, 10, 100, 1000, 10000)) +
  theme_minimal() +
  geom_point(aes(y = n_fire, size = in_interval)) +
  xlab('') +
  ylab('Number of fires > 1000 acres') +
  scale_color_manual(values = c('red', 'black')) +
  scale_size_manual(values = c(1, .3)) +
  theme(legend.position = 'none',
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 6),
        axis.text.x = element_text(angle = 90, size = 8),
        axis.text.y = element_text(size = 8))
ggsave('fig/count-preds.png', width = 10, height = 10)

# overall posterior interval coverage
test_interval_d <- holdout_full_rep %>%
  filter(model == 'zinb_fit.rds') %>%
  left_join(holdout_counts) %>%
  mutate(in_interval = n_fire >= lo & n_fire <= hi)

test_interval_d %>%
  summarize(mean(in_interval))

# worst performance: lowest interval coverage by L3 ecoregion
test_interval_d %>%
  group_by(NA_L3NAME) %>%
  summarize(interval_coverage = mean(in_interval)) %>%
  arrange(interval_coverage)

# what fraction of ecoregions had 100% interval coverage?
test_interval_d %>%
  group_by(NA_L3NAME) %>%
  summarize(interval_coverage = mean(in_interval)) %>%
  left_join(area_df) %>%
  ungroup %>%
  summarize(p_100_pct = mean(interval_coverage == 1),
            n_100_pct = sum(interval_coverage == 1),
            pct_area_100_pct = sum(area[interval_coverage == 1]) / sum(area))

# zoom in on problematic ecoregions
train_full_rep %>%
  filter(model == 'zinb_fit.rds') %>%
  mutate(train = 'train') %>%
  full_join(holdout_full_rep %>% mutate(train = 'test')) %>%
  full_join(er_df) %>%
  filter(NA_L3NAME == 'Cross Timbers') %>%
  left_join(count_df) %>%
  ggplot(aes(ym, med)) +
  geom_point(aes(y = n_fire)) +
  facet_wrap(~NA_L3NAME) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = train), alpha = .4) +
  scale_fill_manual(values = c('red', 'dodgerblue')) +
  #scale_y_log10(breaks = c(0, 1, 10, 100, 1000)) +
  theme_minimal()


