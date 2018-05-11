
library(tidyverse)
library(ggridges)
library(patchwork)
library(sf)

st_covs <- read_rds('data/processed/st_covs.rds')
holdout_counts <- read_rds('data/processed/holdout_counts.rds')

# get areas for each L3 ecoregion
area_df <- read_rds('data/processed/ecoregions.rds') %>%
  as('Spatial') %>%
  as.data.frame %>%
  tbl_df %>%
  group_by(NA_L3NAME) %>%
  summarize(area = sum(Shape_Area))


# Extract posterior draws and visualize results ---------------------------
post <- rstan::extract(read_rds(path = 'zinb_full_fit.rds'),
                       pars = 'count_pred')
str(post)

# turn into data frame
preds <- post$count_pred %>%
  reshape2::melt(varnames = c('iter', 'id')) %>%
  as_tibble %>%
  left_join(select(st_covs, id, ym, NA_L3NAME, year))

# save count predictions as rds
write_rds(preds, 'count-preds.rds')

pred_summary <- preds %>%
  group_by(ym, NA_L3NAME) %>%
  summarize(med = median(value),
            hi = quantile(value, 0.975),
            lo = quantile(value, .025))


# overall posterior interval coverage
test_interval_d <- holdout_counts %>%
  left_join(pred_summary) %>%
  mutate(in_interval = n_fire >= lo & n_fire <= hi)


write_csv(test_interval_d, 'data/processed/count_test_intervals.csv')


# Plot for all ecoregions -------------------------------------------------
test_interval_d %>%
  mutate(in_interval = n_fire >= lo & n_fire <= hi,
         l3_er = gsub(' and ', ' & ', NA_L3NAME),
         longer_than_limit = nchar(l3_er) > 32,
         l3_er = substr(l3_er, start = 1, stop = 32),
         l3_er = ifelse(longer_than_limit,
                        paste0(l3_er, '...'),
                        l3_er)) %>%
  group_by(NA_L3NAME) %>%
  mutate(total_n = sum(n_fire)) %>%
  filter(total_n > 0) %>%
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


#################################################################
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
