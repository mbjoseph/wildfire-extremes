library(tidyverse)
library(rstan)
library(ggridges)
library(viridis)
library(ggrepel)
library(patchwork)
library(ggthemes)

colnamesX <- read_rds('data/processed/colnamesX.rds')
X <- read_rds('data/processed/X.rds')
st_covs <- read_rds('data/processed/st_covs.rds')
ecoregions <- read_rds('data/processed/ecoregions.rds')
stan_d <- read_rds('data/processed/stan_d.rds')

# Extract posterior draws and visualize results ---------------------------
post <- rstan::extract(read_rds('ba_lognormal_fit.rds'), 
                       pars = c('beta', 'lp__', 'mu_full'))
str(post)
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
         nonzero = p_neg > .8 | p_pos > .8,
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

# show important coefficients
beta_summary %>%
  filter(p_neg > .8 | p_pos > .8) %>%
  select(dim, col, p_neg, p_pos, variable,
         median, variable) %>%
  left_join(beta_df) %>%
  ggplot(aes(value, reorder(variable, median), fill = median)) +
  theme_minimal() +
  geom_density_ridges(scale = 3, rel_min_height = 0.005,
                      color = alpha(1, .6)) +
  geom_vline(xintercept = 0, linetype = 'dashed', alpha = .1) +
  theme(axis.text.y = element_text(size = 7)) +
  scale_fill_gradient2() +
  theme(legend.position = 'none',
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(color = 'grey95'),
        axis.text.y = element_text(size = 6)) +
  ylab('') +
  xlab('')

# Partial effect plots ----------------------------------------------------
which_var <- 'rmin'

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
close(pb)

p <- partial_effs %>%
  bind_rows %>%
  left_join(st_covs) %>%
  mutate(NA_L1NAME = tolower(NA_L1NAME),
         NA_L1NAME = factor(tools::toTitleCase(NA_L1NAME)),
         NA_L1NAME = fct_reorder(NA_L1NAME, rmin)) %>%
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

location_ts_plot <- loc_ts %>%
#  filter(year < stan_d$cutoff_year) %>%
  group_by(NA_L1NAME) %>%
  mutate(alpha_val = .5 / length(unique(NA_L3NAME))) %>%
  ggplot(aes(ym, exp(median) + stan_d$min_size,
             fill = NA_L1NAME,
             group = NA_L3NAME)) +
  theme_minimal() +
  geom_ribbon(aes(ymin = exp(lo) + stan_d$min_size,
                  ymax = exp(hi) + stan_d$min_size,
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
  ylab('Expected fire size') +
  facet_wrap(~ NA_L1NAME, nrow = 2,
             labeller = labeller(.rows = label_wrap_gen(23))) +
  ggtitle('C')
location_ts_plot

# plotting expected values vs. covariates
cmap <- c(viridis(6, option = 'C'),
          rev(viridis(6, option = 'C')))

humidity_scatter <- st_covs %>%
  full_join(mu_df) %>%
  filter(year < stan_d$cutoff_year) %>%
  mutate(l2_er = tools::toTitleCase(tolower(as.character(NA_L2NAME))),
         l2_er = gsub(' and ', ' & ', l2_er),
         l2_er = gsub('Usa ', '', l2_er)) %>%
  ggplot(aes(x = rmin,
             y = exp(median) + stan_d$min_size,
             color = month)) +
  geom_point(alpha = .6) +
  theme_minimal() +
  scale_color_gradientn(colors = cmap, 'Month') +
  xlab('Mean daily minimum humidity') +
  ylab('Expected fire size') +
  theme(panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 8, color = 'grey30')) +
  ggtitle('B')
humidity_scatter

q <- (coefplot + p) / location_ts_plot + plot_layout(heights = c(.6, 1), ncol = 1)
q
ggsave('fig/burn-area-effs.png', plot = q, width = 8.5, height = 5)
