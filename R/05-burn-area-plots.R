source('R/02-explore.R')
soucr('R/make-stan-d.R')
library(extraDistr)
library(ggridges)
library(ggthemes)
library(plotly)
library(viridis)
library(ggrepel)

fit <- read_rds(path = list.files(pattern = 'lognormal_.*'))

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
beta_summary %>%
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
                  alpha = 1, size = 2.5)

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


beta_summary %>%
  filter(p_neg > .8 | p_pos > .8)

gc()

# Visualize some predictions ----------------------------------------------
st_covs$row <- 1:nrow(st_covs)

mu_df <- post$mu_full %>%
  reshape2::melt(varnames = c('iter', 'row')) %>%
  tbl_df %>%
  group_by(row) %>%
  summarize(median = median(value),
            lo = quantile(value, 0.05),
            hi = quantile(value, 0.95),
            q1 = quantile(value, .1),
            q2 = quantile(value, .2),
            q3 = quantile(value, .3),
            q4 = quantile(value, .4),
            q6 = quantile(value, .6),
            q7 = quantile(value, .7),
            q8 = quantile(value, .8),
            q9 = quantile(value, .9)) %>%
  ungroup

## Visualize parameter time series
plot_mu_ts <- function(df) {
  df %>%
    ggplot(aes(ym, median, fill = NA_L1NAME)) +
    theme_minimal() +
    geom_ribbon(aes(ymin = lo, ymax = hi),
                alpha = .3, color = NA) +
    geom_ribbon(aes(ymin = q1, ymax = q9),
                alpha = .3, color = NA) +
    geom_ribbon(aes(ymin = q2, ymax = q8),
                alpha = .3, color = NA) +
    geom_ribbon(aes(ymin = q3, ymax = q7),
                alpha = .3, color = NA) +
    geom_ribbon(aes(ymin = q4, ymax = q6),
                alpha = .3, color = NA) +
    geom_line(size = .1) +
    facet_wrap(~ facet_factor) +
    xlab('Date') +
    geom_vline(xintercept = cutoff_year,
               linetype = 'dashed', col = 'grey') +
    scale_color_gdocs() +
    scale_fill_gdocs() +
    theme(legend.position = 'none')
}

# location parameter
mu_df %>%
  full_join(st_covs) %>%
  left_join(dplyr::distinct(tbl_df(ecoregions),
                            NA_L3NAME, NA_L1CODE)) %>%
  mutate(facet_factor = paste(NA_L1CODE, NA_L3NAME)) %>%
  plot_mu_ts +
  ylab('Lognormal mean')

# plotting expected values vs. covariates
cmap <- c(viridis(12, option = 'C'),
          rev(viridis(12, option = 'C')))

mu_df %>%
  full_join(st_covs) %>%
  left_join(dplyr::distinct(tbl_df(ecoregions),
                            NA_L3NAME, NA_L1CODE)) %>%
  mutate(facet_factor = paste(NA_L1CODE, NA_L3NAME)) %>%
  ggplot(aes(rmin, exp(median), color = month)) +
  scale_color_gradientn(colors = cmap) +
  theme_minimal() +
  facet_wrap(~NA_L2NAME) +
  xlab('Mean minimum daily relative humidity') +
  ylab('Expected burn area') +
  scale_y_log10() +
  geom_linerange(aes(ymin = exp(q1), ymax = exp(q9)), alpha = .5) +
  geom_point(alpha = .9, size = .5)



st_covs %>%
  full_join(mu_df) %>%
  filter(year < cutoff_year) %>%
  mutate(l2_er = tools::toTitleCase(tolower(as.character(NA_L2NAME))),
         l2_er = gsub(' and ', ' & ', l2_er),
         l2_er = gsub('Usa ', '', l2_er)) %>%
  ggplot(aes(x = rmin,
             y = exp(median),
             color = month)) +
  geom_linerange(aes(ymin = exp(q1),
                     ymax = exp(q9)),
                 alpha = .1) +
  geom_point(size = .5) +
  theme_minimal() +
  facet_wrap(~ fct_reorder(l2_er, rmin),
             labeller = labeller(.rows = label_wrap_gen(25))) +
  scale_color_gradientn(colors = cmap, 'Month') +
  scale_y_log10() +
  xlab('Mean daily minimum humidity') +
  ylab('Expected burn area exceedance (acres > 1000)') +
  theme(panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 8, color = 'grey30'))
ggsave('fig/humidity-burn-area.png', width = 9, height = 6)
