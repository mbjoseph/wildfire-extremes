library(tidyverse)
library(rstan)
library(viridis)
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
  select(-year)

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
  as_tibble

ln_sigma <- ln_post$scale[1:count_post_iters] %>%
  reshape2::melt(varnames = c('iter'), value.name = 'ln_scale') %>%
  as_tibble

# merge lognormal mean and scale draws with counts
test_preds <- full_join(zinb_preds, ln_mu) %>%
  left_join(ln_sigma)

# for every event predicted by the count model, simulate burn areas
test_preds <- test_preds %>%
  filter(n_event > 0) %>%
  rowwise() %>%
  mutate(preds = list(exp(rnorm(n_event, ln_mu, ln_scale))),
         pred_max = max(unlist(preds)),
         tot_max = sum(unlist(preds)))

test_pred_summary <- test_preds %>%
  ungroup %>%
  group_by(NA_L3NAME, ym) %>%
  summarize(med_max = median(pred_max),
            lo_max = quantile(pred_max, .025),
            hi_max = quantile(pred_max, .975),
            med_sum = median(tot_max),
            lo_sum = quantile(tot_max, .025),
            hi_sum = quantile(tot_max, .975)) %>%
  ungroup()


cmap <- c(viridis(12, option = 'C'),
          rev(viridis(12, option = 'C')))


test_pred_summary %>%
  left_join(select(st_covs, NA_L3NAME, ym, rmin,
                   month, NA_L2NAME)) %>%
  ggplot(aes(rmin, med_max, color = month)) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank()) +
  geom_linerange(aes(ymin = lo_max, ymax = hi_max), alpha = .05) +
  geom_point() +
  scale_y_log10() +
  scale_color_gradientn(colors = cmap, 'Month') +
  xlab('Mean daily minimum relative humidity') +
  ylab('Predicted maximum burn area exceedance')

test_pred_summary %>%
  left_join(select(st_covs, NA_L3NAME, ym, rmin,
                   month, NA_L2NAME)) %>%
  ggplot(aes(rmin, med_sum, color = month)) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank()) +
  geom_linerange(aes(ymin = lo_sum, ymax = hi_sum), alpha = .05) +
  geom_point() +
  scale_y_log10() +
  scale_color_gradientn(colors = cmap, 'Month') +
  xlab('Mean daily minimum relative humidity') +
  ylab('Predicted total burn area exceedance')

test_pred_summary %>%
  ggplot(aes(x = med_sum, y = med_max)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10()
