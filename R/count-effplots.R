
# Count plots -------------------------------------------------------------

source('R/02-explore.R')
library(extraDistr)
library(ggridges)
library(viridis)
library(ggthemes)

fit <- read_rds(path = list.files(pattern = 'nb_fit.*'))

# Evaluate convergence ----------------------------------------------------
traceplot(fit, inc_warmup = TRUE)

traceplot(fit, pars = c('tau', 'c', 'alpha'))

# Extract posterior draws and visualize results ---------------------------
post <- rstan::extract(fit)
str(post)
rm(fit)
gc()



# Coefficients that seem to be far from zero ------------------------------
beta_df <- post$beta %>%
  reshape2::melt(varnames = c('iter', 'dim', 'col')) %>%
  tbl_df

beta_summary <- beta_df %>%
  group_by(dim, col) %>%
  summarize(median = median(value),
            lo = quantile(value, 0.025),
            hi = quantile(value, 0.975),
            p_neg = mean(value < 0),
            p_pos = mean(value > 0)) %>%
  ungroup %>%
  mutate(variable = colnamesX[col],
         nonzero = p_neg > .9 | p_pos > .9)

beta_summary %>%
  filter(nonzero) %>%
  select(dim, col, p_neg, p_pos, variable) %>%
  left_join(beta_df) %>%
  ggplot(aes(value, variable)) +
  geom_density_ridges() +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  facet_wrap(~ dim, scales = 'free') +
  theme(axis.text.y = element_text(size = 7)) +
  ggtitle('Estimated effects on # of fires')
ggsave('fig/fire-effs.png', width = 8, height = 8)

beta_summary %>%
  filter(nonzero) %>%
  ggplot(aes(x = median, y = reorder(variable, median))) +
  geom_point() +
  geom_errorbarh(aes(xmin = lo, xmax = hi)) +
  theme(axis.text.y = element_text(size = 5))
ggsave('fig/count-effs.png', width = 6, height = 12)

rm(beta_df)
gc()

mu_df <- post$mu_full %>%
  reshape2::melt(varnames = c('iter', 'dim', 'id')) %>%
  as_tibble() %>%
  group_by(id) %>%
  summarize(med = median(value),
            lo = quantile(value, .05),
            hi = quantile(value, .95))

st_covs <- st_covs %>%
  group_by(NA_L3NAME) %>%
  mutate(mean_hd = mean(housing_density))


st_covs %>%
  full_join(mu_df) %>%
  ggplot(aes(x = rmin, y = exp(med), color = housing_density)) +
  geom_point(size = .5, alpha = .5) +
  theme_minimal() +
  facet_wrap(~ reorder(NA_L3NAME, mean_hd)) +
  scale_color_viridis(trans = 'log') +
  scale_y_log10()

st_covs %>%
  full_join(mu_df) %>%
  ggplot(aes(x = vs, y = med, color = housing_density)) +
  geom_point(size = .5, alpha = .5) +
  theme_minimal() +
  facet_wrap(~ reorder(NA_L3NAME, mean_hd)) +
  scale_color_viridis(trans = 'log')

st_covs %>%
  full_join(mu_df) %>%
  ggplot(aes(x = prev_12mo_precip, y = med, color = housing_density)) +
  geom_point(size = .5, alpha = .5) +
  theme_minimal() +
  facet_wrap(~ reorder(NA_L3NAME, mean_hd)) +
  scale_color_viridis(trans = 'log')

st_covs %>%
  full_join(mu_df) %>%
  ggplot(aes(x = pr, y = med, color = housing_density)) +
  geom_point(size = .5, alpha = .5) +
  theme_minimal() +
  facet_wrap(~ reorder(NA_L3NAME, mean_hd)) +
  scale_color_viridis(trans = 'log')

st_covs %>%
  full_join(mu_df) %>%
  ggplot(aes(x = tmmx, y = exp(med), color = housing_density)) +
  geom_point(size = .5, alpha = .5) +
  theme_minimal() +
  facet_wrap(~ reorder(NA_L3NAME, mean_hd)) +
  scale_color_viridis(trans = 'log') +
  scale_y_log10()
