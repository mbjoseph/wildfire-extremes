library(tidyverse)
library(rstan)
library(ggrepel)
library(ggthemes)
library(patchwork)

colnamesX <- read_rds('data/processed/colnamesX.rds')
X <- read_rds('data/processed/X.rds')
st_covs <- read_rds('data/processed/st_covs.rds')
ecoregions <- read_rds('data/processed/ecoregions.rds')
stan_d <- read_rds('data/processed/stan_d.rds')

# Extract posterior draws and visualize results ---------------------------
post <- rstan::extract(read_rds('ba_lognormal_fit.rds'), 
                       pars = c('beta', 'lp__', 'mu_full'))
str(post)

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
         nonzero = p_neg > .8 | p_pos > .8,
         variable = gsub('bs_', '', variable),
         variable = gsub('rmin', 'humidity', variable),
         variable = gsub('tmmx', 'temperature', variable),
         variable = gsub('vs', 'wind speed', variable),
         variable = gsub('prev_12mo_precip', '12 mo. precip.', variable),
         variable = gsub('pr', 'precipitation', variable),
         varirable = gsub('log_housing_density', 'housing density', variable),
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
coefplot <- beta_summary %>%
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
  geom_segment(aes(y = lo, yend = hi,
                   xend = variable),
               size = .3) +
  xlab('Coefficient') +
  ylab('Coefficient value') +
  scale_color_gradient2(mid = 'black', low = 'blue', high = 'red') +
  theme(legend.position = 'none') +
  geom_text_repel(aes(label = ifelse(abs(rel_size) > .35, variable, '')),
                  alpha = 1, size = 2.5) +
  ggtitle('A')
coefplot

write_csv(x = beta_summary, path = 'data/processed/burn-area-beta.csv')

# Partial effect plots ----------------------------------------------------
which_var <- c('rmin', 'tmmx', 'vs', 'prev_12mo_precip', 'log_housing_density')

partial_effs_each <- vector(mode = 'list', length = length(which_var))
for (k in seq_along(which_var)) {
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
    cols <- grepl(which_var[k], colnames(X_sub))
    
    effects <- array(dim = c(nrow(X_sub), 3)) # 3: med, lo, hi
    for (j in 1:nrow(X_sub)) {  # month j
      vals <- X_sub[j, cols] %*% t(post$beta[, 1, cols])
      effects[j, 1] <- quantile(vals, .025)
      effects[j, 2] <- median(vals)
      effects[j, 3] <- quantile(vals, .975)
    }
    partial_effs[[i]] <- effects %>%
      reshape2::melt(varnames = c('row_id', 'quantity')) %>%
      as_tibble %>%
      mutate(quantity = case_when(.$quantity == 1 ~ 'lo',
                                  .$quantity == 2 ~ 'med',
                                  .$quantity == 3 ~ 'hi')) %>%
      spread(quantity, value) %>%
      left_join(select(subdf, row_id, NA_L3NAME, year, ym, housing_density))
  }
  partial_effs_each[[k]] <- partial_effs
  close(pb)
}

names(partial_effs_each) <- which_var

write_rds(partial_effs_each, 'data/processed/burn_area_partial_all.rds')


burn_area_partials <- partial_effs_each[['rmin']] %>%
  bind_rows %>%
  left_join(st_covs) %>%
  mutate(NA_L1NAME = tolower(NA_L1NAME),
         NA_L1NAME = factor(tools::toTitleCase(NA_L1NAME)),
         NA_L1NAME = fct_reorder(NA_L1NAME, rmin))

write_rds(burn_area_partials, 'data/processed/burn_area_partials.rds')

p <- burn_area_partials %>%
  ggplot(aes(rmin, med, group = NA_L3NAME)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), color = NA,
              alpha = .1) +
  scale_fill_gdocs() +
  geom_line(alpha = .3) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        legend.position = 'none') +
  xlab('Mean daily minimum humidity') +
  ylab('Partial effect') +
  ggtitle('B')
p




# Visualize some predictions ----------------------------------------------
st_covs$row <- 1:nrow(st_covs)

mu_df <- post$mu_full %>%
  reshape2::melt(varnames = c('iter', 'row')) %>%
  tbl_df %>%
  group_by(row) %>%
  summarize(median = median(value),
            lo = quantile(value, 0.05),
            hi = quantile(value, 0.95)) %>%
  ungroup

# location parameter
loc_ts <- mu_df %>%
  full_join(st_covs) %>%
  left_join(dplyr::distinct(tbl_df(ecoregions),
                            NA_L3NAME, NA_L1CODE)) %>%
  mutate(facet_factor = paste(NA_L1CODE, NA_L3NAME)) %>%
  mutate(NA_L1NAME = tolower(NA_L1NAME),
         NA_L1NAME = factor(tools::toTitleCase(NA_L1NAME)),
         NA_L1NAME = fct_reorder(NA_L1NAME, rmin))

write_rds(loc_ts, 'data/processed/loc_ts.rds')

hectares_per_acre <- 0.404686
location_ts_plot <- loc_ts %>%
  group_by(NA_L1NAME) %>%
  mutate(alpha_val = .5 / length(unique(NA_L3NAME))) %>%
  ggplot(aes(ym, (exp(median) + stan_d$min_size) * hectares_per_acre,
             fill = NA_L1NAME,
             group = NA_L3NAME)) +
  theme_minimal() +
  geom_ribbon(aes(ymin = (exp(lo) + stan_d$min_size) * hectares_per_acre,
                  ymax = (exp(hi) + stan_d$min_size) * hectares_per_acre,
                  alpha = alpha_val),
              color = NA) +
  geom_line(size = .2, aes(alpha = alpha_val)) +
  xlab('') +
  scale_alpha_continuous(range = c(.2, .6)) +
  scale_color_gdocs() +
  scale_fill_gdocs() +
  theme(legend.position = 'none',
        panel.grid.minor = element_blank()) +
  scale_y_log10() +
  ylab('Expected fire size (hectares)') +
  facet_wrap(~ NA_L1NAME, nrow = 2,
             labeller = labeller(.rows = label_wrap_gen(23))) +
  ggtitle('C')
location_ts_plot

q <- (coefplot + p) / location_ts_plot + plot_layout(heights = c(.6, 1), ncol = 1)
q
ggsave('fig/figure_8.pdf', plot = q, width = 8.5, height = 5)
