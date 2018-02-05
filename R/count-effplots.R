
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
            mean = mean(value),
            lo = quantile(value, 0.025),
            hi = quantile(value, 0.975),
            p_neg = mean(value < 0),
            p_pos = mean(value > 0)) %>%
  ungroup %>%
  mutate(variable = colnamesX[col],
         nonzero = p_neg > .95 | p_pos > .95)

beta_summary %>%
  filter(nonzero) %>%
  select(dim, col, p_neg, p_pos, variable, median, mean) %>%
  left_join(beta_df) %>%
  mutate(Response = ifelse(dim == 1,
                           'Negative binomial mean',
                           'Zero-inflation component'),
         variable = gsub('bs_', '', variable),
         variable = gsub('crmin', 'humidity', variable),
         variable = gsub('ctmx', 'temperature', variable),
         variable = gsub('cvs', 'wind speed', variable),
         variable = gsub('cpr12', '12 mo. precip.', variable),
         variable = gsub('cpr', 'precipitation', variable),
         variable = ifelse(grepl(':', x = variable),
                           paste0('Intxn(', variable, ')'),
                           variable),
         variable = gsub(':', ' x ', variable),
         variable = gsub('NA_', '', variable),
         variable = gsub('NAME', ' ', variable),
         variable = gsub('_', ': ', variable),
         variable = tolower(variable),
         variable = tools::toTitleCase(variable)) %>%
  ggplot(aes(value, reorder(variable, mean), fill = median)) +
  theme_minimal() +
  geom_density_ridges(scale = 3, rel_min_height = 0.005,
                      color = alpha(1, .6)) +
  geom_vline(xintercept = 0, linetype = 'dashed', alpha = .1) +
  facet_wrap(~ Response, scales = 'free', nrow = 2, shrink = TRUE) +
  theme(axis.text.y = element_text(size = 7)) +
  scale_fill_gradient2() +
  theme(legend.position = 'none',
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(color = 'grey95'),
        axis.text.y = element_text(size = 6)) +
  ylab('') +
  xlab('')
ggsave('fig/fire-effs.pdf', width = 6, height = 7)

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
         facet_factor = paste(NA_L2CODE, NA_L3NAME, sep = ': '),
         l2_er = tools::toTitleCase(tolower(as.character(NA_L2NAME))),
         l2_er = gsub(' and ', ' & ', l2_er),
         l2_er = gsub('Usa ', '', l2_er))

cmap <- c(viridis(12, option = 'C'),
          rev(viridis(12, option = 'C')))

st_covs %>%
  full_join(mu_df) %>%
  filter(year < cutoff_year) %>%
  ggplot(aes(x = rmin,
             y = med / (1000 * area),
             color = month)) +
  geom_linerange(aes(ymin = lo / (1000 * area),
                     ymax = hi / (1000 * area)),
                 alpha = .1) +
  geom_point(size = .5) +
  theme_minimal() +
  facet_wrap(~ fct_reorder(l2_er, rmin),
             labeller = labeller(.rows = label_wrap_gen(25))) +
  scale_color_gradientn(colors = cmap, 'Month') +
  scale_y_log10() +
  xlab('Mean daily minimum humidity') +
  ylab(expression(paste('Expected fire density: # per ', km^2))) +
  theme(panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 8, color = 'grey30'))
ggsave('fig/humidity-counts.png', width = 9, height = 6)

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

