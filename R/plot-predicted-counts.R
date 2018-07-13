
library(tidyverse)
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
write_csv(area_df, 'data/processed/area_df.csv')

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
