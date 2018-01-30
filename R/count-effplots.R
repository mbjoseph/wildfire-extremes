
# Count plots -------------------------------------------------------------

source('R/02-explore.R')
source('R/make-stan-d.R')
library(ggridges)
library(viridis)
library(ggthemes)
library(plotly)

fit <- read_rds(path = list.files(pattern = 'zinb_fit.*'))

# Evaluate convergence ----------------------------------------------------
traceplot(fit, inc_warmup = TRUE)

traceplot(fit, pars = c('tau', 'c', 'alpha', 'Rho_beta'))

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
  select(dim, col, p_neg, p_pos, variable, median) %>%
  left_join(beta_df) %>%
  mutate(Response = ifelse(dim == 1,
                           'Neg. binom mean',
                           'Possibility of fire')) %>%
  ggplot(aes(value, reorder(variable, median), fill = median)) +
  theme_classic() +
  geom_density_ridges(scale = 10, color = alpha(1, .1)) +
  geom_vline(xintercept = 0, linetype = 'dashed', alpha = .1) +
  facet_wrap(~ Response, scales = 'free') +
  theme(axis.text.y = element_text(size = 7)) +
  ggtitle('Estimated effects on # of fires') +
  scale_fill_gradient2() +
  theme(legend.position = 'none')
ggsave('fig/fire-effs.png', width = 12, height = 8)

beta_summary %>%
  filter(nonzero) %>%
  ggplot(aes(x = median, y = reorder(variable, median))) +
  geom_point() +
  facet_wrap(~dim, scales = 'free') +
  geom_errorbarh(aes(xmin = lo, xmax = hi)) +
  theme(axis.text.y = element_text(size = 8))
ggsave('fig/count-effs.png', width = 6, height = 12)

rm(beta_df)
gc()

mu_df <- post$mu_full %>%
  reshape2::melt(varnames = c('iter', 'dim', 'id')) %>%
  as_tibble() %>%
  spread(dim, value) %>%
  mutate(expected_value = exp(`1`) * plogis(`2`)) %>%
  group_by(id) %>%
  summarize(med = median(expected_value),
            lo = quantile(expected_value, .05),
            hi = quantile(expected_value, .95))

st_covs <- st_covs %>%
  group_by(NA_L3NAME) %>%
  mutate(mean_hd = mean(housing_density),
         facet_factor = paste(NA_L2CODE, NA_L3NAME, sep = ': '))

cmap <- c(viridis(12, option = 'C'),
          rev(viridis(12, option = 'C')))

st_covs %>%
  full_join(mu_df) %>%
  filter(year < cutoff_year) %>%
  ggplot(aes(x = rmin,
             y = med / (1000 * area),
             color = month)) +
  geom_point(size = .5) +
  theme_minimal() +
  facet_wrap(~ NA_L2NAME,
             labeller = labeller(NA_L2NAME = label_wrap_gen(20))) +
  scale_color_gradientn(colors = cmap) +
  scale_y_log10() +
  xlab('Mean daily minimum humidity') +
  ylab('Expected fire density (# per sq. km)')

st_covs %>%
  full_join(mu_df) %>%
  filter(year < cutoff_year) %>%
  ggplot(aes(x = vs,
             y = med / (1000 * area),
             color = month)) +
  geom_point(size = .5) +
  theme_minimal() +
  facet_wrap(~ NA_L2NAME,
             labeller = labeller(NA_L2NAME = label_wrap_gen(20))) +
  scale_color_gradientn(colors = cmap) +
  scale_y_log10() +
  xlab('Mean daily wind speed') +
  ylab('Expected fire density (# per sq. km)')

st_covs %>%
  full_join(mu_df) %>%
  filter(year < cutoff_year) %>%
  ggplot(aes(x = prev_12mo_precip,
             y = med / (1000 * area),
             color = month)) +
  geom_point(size = .5) +
  theme_minimal() +
  facet_wrap(~ NA_L2NAME,
             labeller = labeller(NA_L2NAME = label_wrap_gen(20))) +
  scale_color_gradientn(colors = cmap) +
  scale_y_log10() +
  xlab('Mean previous 12 month precipitation') +
  ylab('Expected fire density (# per sq. km)') +
  scale_x_log10()

st_covs %>%
  full_join(mu_df) %>%
  filter(year < cutoff_year) %>%
  ggplot(aes(x = pr,
             y = med / (1000 * area),
             color = month)) +
  geom_point(size = .5) +
  theme_minimal() +
  facet_wrap(~ NA_L2NAME,
             labeller = labeller(NA_L2NAME = label_wrap_gen(20))) +
  scale_color_gradientn(colors = cmap) +
  scale_y_log10() +
  xlab('Same month precipitation') +
  ylab('Expected fire density (# per sq. km)') +
  scale_x_log10()

st_covs %>%
  full_join(mu_df) %>%
  filter(year < cutoff_year) %>%
  ggplot(aes(x = tmmx,
             y = med / (1000 * area),
             color = month,
             group = NA_L3NAME)) +
  geom_point(size = .5) +
  theme_minimal() +
  facet_wrap(~ NA_L2NAME,
             labeller = labeller(NA_L2NAME = label_wrap_gen(20))) +
  scale_color_gradientn(colors = cmap) +
  scale_y_log10() +
  xlab('Air temperature') +
  ylab('Expected fire density (# per sq. km)')

