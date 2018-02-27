
# Count plots -------------------------------------------------------------

source('R/02-explore.R')
source('R/make-stan-d.R')
library(ggridges)
library(viridis)
library(ggthemes)
library(plotly)
library(ggrepel)

fit <- read_rds(path = 'zinb_full_fit.rds')

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
         nonzero = p_neg > .85 | p_pos > .85) %>%
  mutate(Response = ifelse(dim == 1,
                           'Negative binomial component',
                           'Zero-inflation component'),
         variable = gsub('bs_', '', variable),
         variable = gsub('crmin', 'humidity', variable),
         variable = gsub('ctmx', 'temperature', variable),
         variable = gsub('cvs', 'wind speed', variable),
         variable = gsub('cpr12', '12 mo. precip.', variable),
         variable = gsub('cpr', 'precipitation', variable),
         variable = gsub('chd', 'Housing density', variable),
         variable = ifelse(grepl(':', x = variable),
                           paste0('Intxn(', variable, ')'),
                           variable),
         variable = gsub(':', ' x ', variable),
         variable = gsub('NA_', '', variable),
         variable = gsub('NAME', ' ', variable),
         variable = gsub('_', ': ', variable),
         variable = tolower(variable),
         variable = tools::toTitleCase(variable))

# show all coefficients
beta_summary %>%
  group_by(Response) %>%
  mutate(max_abs = max(abs(median)),
         rel_size = median / max_abs) %>%
  ggplot(aes(y = median,
             x = variable,
             color = rel_size ^ 2 * sign(rel_size),
             alpha = rel_size ^ 2)) +
  geom_point(size = .3) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank()) +
  facet_wrap(~ Response, scales = 'free', nrow = 2) +
  geom_segment(aes(y = lo, yend = hi,
                   xend = variable),
               size = .3) +
  xlab('Coefficient') +
  ylab('Coefficient value') +
  scale_color_gradient2(mid = 'black', low = 'blue', high = 'red') +
  theme(legend.position = 'none') +
  geom_text_repel(aes(label = ifelse(abs(rel_size) > .35, variable, '')),
                  alpha = 1, size = 2.5)
ggsave('fig/all-coefs.png', width = 8, height = 4)


# show important coefficients
beta_summary %>%
  filter(nonzero) %>%
  select(dim, col, p_neg, p_pos, variable,
         median, mean, variable, Response) %>%
  left_join(beta_df) %>%
  ggplot(aes(value, reorder(variable, mean), fill = median)) +
  theme_minimal() +
  geom_density_ridges(scale = 3, rel_min_height = 0.005,
                      color = alpha(1, .6)) +
  geom_vline(xintercept = 0, linetype = 'dashed', alpha = .1) +
  facet_wrap(~ Response, scales = 'free',
             shrink = TRUE, ncol = 1) +
  theme(axis.text.y = element_text(size = 7)) +
  scale_fill_gradient2() +
  theme(legend.position = 'none',
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(color = 'grey95')) +
  ylab('') +
  xlab('')
ggsave('fig/fire-effs.pdf', width = 6, height = 7)

rm(beta_df)
gc()


# Partial effect plots ----------------------------------------------------
which_var <- c('ctmx')

partial_effs <- list()
n_iter <- length(post$lp__)
unique_ers <- unique(st_covs$NA_L3NAME)
pb <- txtProgressBar(max = length(unique_ers), style = 3)
for (i in seq_along(unique_ers)) {
  setTxtProgressBar(pb, i)
  subdf <- st_covs %>%
    filter(NA_L3NAME == unique_ers[i]) %>%
    mutate(row_id = 1:n())
  X_sub <- X[st_covs$NA_L3NAME == unique_ers[i], ]
  cols <- grepl(which_var, colnames(X_sub))

  effects <- array(dim = c(nrow(X_sub), 2, 3)) # 2 responses, 3: med, lo, hi
  for (j in 1:nrow(X_sub)) {  # month j
    for (k in 1:2) {          # response k
      vals <- X_sub[j, cols] %*% t(post$beta[1:100, k, cols])
      effects[j, k, 1] <- quantile(vals, .025)
      effects[j, k, 2] <- median(vals)
      effects[j, k, 3] <- quantile(vals, .975)
    }
  }
  partial_effs[[i]] <- effects %>%
    reshape2::melt(varnames = c('row_id', 'response', 'quantity')) %>%
    as_tibble %>%
    mutate(response = ifelse(response == 1, 'negbinom', 'zeroinfl'),
           quantity = case_when(.$quantity == 1 ~ 'lo',
                                .$quantity == 2 ~ 'med',
                                .$quantity == 3 ~ 'hi')) %>%
    spread(quantity, value) %>%
    left_join(select(subdf, row_id, NA_L3NAME, year, ym, housing_density))
}
close(pb)

p <- partial_effs %>%
  bind_rows %>%
  left_join(st_covs) %>%
  filter(response == 'negbinom') %>%
  ggplot(aes(tmmx, med, group = NA_L3NAME, color = NA_L1NAME)) +
#  geom_ribbon(aes(ymin = lo, ymax = hi), color = NA, alpha = .1) +
  geom_line(alpha = .8) +
  theme_minimal() +
  #facet_wrap(~ NA_L1NAME, scales = 'free_y') +
  #scale_x_log10() +
  xlab('Mean daily temperature') +
  ylab('Partial effect')
ggplotly(p)

















# Expected values vs. covariate values ------------------------------------
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
  ylab(expression(paste('Expected fire density: # per (month x ', km^2, ')'))) +
  theme(panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 8, color = 'grey30'))
ggsave('fig/humidity-counts.png', width = 9, height = 6)


# Housing density?? -------------------------------------------------------
p <- st_covs %>%
  full_join(mu_df) %>%
  filter(year < cutoff_year) %>%
  group_by(NA_L3NAME, year) %>%
  mutate(housing_density = median(housing_density),
         med = median(med)) %>%
  distinct(housing_density, med, NA_L3NAME, area) %>%
  ggplot(aes(x = housing_density,
             y = med / (1000 * area),
             group = NA_L3NAME)) +
  geom_point(size = .5, alpha = .4) +
  theme_minimal() +
  scale_y_log10() +
  xlab('Housing density') +
  ylab('Median annual fire density') +
  theme(panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 8, color = 'grey30')) +
  geom_smooth(method = 'lm', se = FALSE) +
  scale_x_log10()
ggplotly(p)


st_covs %>%
  full_join(mu_df) %>%
  filter(year < cutoff_year) %>%
  mutate(temperature_c = tmmx - 273.15) %>%
  ggplot(aes(x = temperature_c,
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
  xlab('Mean daily maximum air temperature (C)') +
  ylab(expression(paste('Expected fire density: # per (month x ', km^2, ')'))) +
  theme(panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 8, color = 'grey30'))
ggsave('fig/air-temp-counts.png', width = 9, height = 6)



st_covs %>%
  full_join(mu_df) %>%
  filter(year < cutoff_year) %>%
  ggplot(aes(x = vs,
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
  xlab('Mean monthly wind speed') +
  ylab(expression(paste('Expected fire density: # per (month x ', km^2, ')'))) +
  theme(panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 8, color = 'grey30'))
ggsave('fig/wind-counts.png', width = 9, height = 6)


st_covs %>%
  full_join(mu_df) %>%
  filter(year < cutoff_year) %>%
  ggplot(aes(x = pr,
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
  xlab('Mean monthly precipitation') +
  ylab(expression(paste('Expected fire density: # per (month x ', km^2, ')'))) +
  theme(panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 8, color = 'grey30')) +
  scale_x_log10()
ggsave('fig/precip-counts.png', width = 9, height = 6)


st_covs %>%
  full_join(mu_df) %>%
  filter(year < cutoff_year) %>%
  ggplot(aes(x = prev_12mo_precip,
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
  xlab('Mean precipitation over previous 12 months') +
  ylab(expression(paste('Expected fire density: # per (month x ', km^2, ')'))) +
  theme(panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 8, color = 'grey30')) +
  scale_x_log10()
ggsave('fig/precip12-counts.png', width = 9, height = 6)

