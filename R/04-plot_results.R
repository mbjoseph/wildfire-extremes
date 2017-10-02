library(raster)
library(rstan)
library(sf)
library(rmapshaper)
library(tidyverse)
library(rstan)
library(modeest)
library(HDInterval)

source('R/02-explore.R')


m_fit <- read_rds('m_fit.rds')




# Evaluate convergence ----------------------------------------------------
traceplot(m_fit, inc_warmup = TRUE)

traceplot(m_fit, pars = c('tau', 'sigma', 'c_sq', 'alpha'))

plot(m_fit, pars = 'beta') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  coord_flip()


# Extract posterior draws and visualize results ---------------------------
post <- rstan::extract(m_fit)
str(post)
rm(m_fit)
gc()



# Coefficients that seem to be far from zero ------------------------------
beta_df <- post$beta %>%
  reshape2::melt(varnames = c('iter', 'dim', 'col')) %>%
  tbl_df %>%
  group_by(dim, col) %>%
  summarize(mode = mlv(value, method = 'venter')$M,
            lo = hdi(value, credMass = .9)[1],
            hi = hdi(value, credMass = .9)[2]) %>%
  ungroup %>%
  mutate(variable = colnames(X)[col],
         over_zero = lo < 0 & hi > 0)

beta_df %>%
  filter(!over_zero) %>%
  ggplot(aes(mode, variable)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_segment(aes(x = lo, xend = hi, yend = variable)) +
  facet_wrap(~ dim, scales = 'free')

beta_df %>%
  ggplot(aes(mode, variable)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_segment(aes(x = lo, xend = hi, yend = variable)) +
  facet_wrap(~ dim, scales = 'free_x')



# explore the possibility of allowing correlation between beta and beta_c
beta_df %>%
  select(-lo, -hi, -over_zero) %>%
  spread(dim, mode) %>%
  ggplot(aes(`1`, `3`)) +
  geom_point(alpha = .5) +
  geom_vline(xintercept = 0,linetype = 'dashed') +
  geom_hline(yintercept = 0,linetype = 'dashed') +
  xlab('Effect on number of events') +
  ylab('Effect on event size')
# does seem to be evidence for correpsondence - maybe multivariate horseshoe

# Visualize some predictions ----------------------------------------------

burn_covs$row <- 1:nrow(burn_covs)

mu_df <- post$mu %>%
  reshape2::melt(varnames = c('iter', 'response', 'row')) %>%
  tbl_df %>%
  group_by(response, row) %>%
  summarize(mode = mlv(value, method = 'venter')$M,
            lo = hdi(value, credMass = .9)[1],
            hi = hdi(value, credMass = .9)[2]) %>%
  ungroup

## Visualize expected burn area exceedance values
burn_mu_df <- mu_df %>%
  filter(response == 1) %>%
  full_join(tbl_df(burn_covs)) %>%
  left_join(dplyr::distinct(tbl_df(ecoregions),
                            NA_L3NAME, NA_L1CODE)) %>%
  mutate(facet_factor = paste0(NA_L1CODE, NA_L3NAME))

# changes over time
plot_mu_ts <- function(df) {
  df %>%
    ggplot(aes(ym, exp(mode), fill = NA_L1NAME)) +
    theme_minimal() +
    geom_ribbon(aes(ymin = exp(lo), ymax = exp(hi)),
                alpha = .5, color = NA) +
    geom_line(size = .1) +
    facet_wrap(~ NA_L3NAME) +
    xlab('Date') +
    ylab("Expected burn area log exceedance") +
    #scale_y_log10() +
    geom_vline(xintercept = cutoff_year, linetype = 'dashed', col = 'grey') +
    scale_color_gdocs() +
    scale_fill_gdocs() +
    theme(legend.position = 'none')
}


burn_mu_df %>%
  filter(NA_L1NAME == 'EASTERN TEMPERATE FORESTS') %>%
  plot_mu_ts

burn_mu_df %>%
  filter(NA_L1NAME == 'NORTH AMERICAN DESERTS') %>%
  plot_mu_ts

burn_mu_df %>%
  filter(NA_L1NAME == 'GREAT PLAINS') %>%
  plot_mu_ts

burn_mu_df %>%
  filter(NA_L1NAME == 'NORTHWESTERN FORESTED MOUNTAINS') %>%
  plot_mu_ts




# true values vs. predicted means
burn_mu_df %>%
  right_join(tbl_df(train_burns)) %>%
  ggplot(aes(log(R_ACRES - 1e3), mode)) +
  geom_point(size = .2, alpha = .2) +
  geom_segment(aes(xend = log(R_ACRES - 1e3),
                   y = lo, yend = hi), alpha = .2) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
  facet_wrap(~ NA_L3NAME) +
  xlab('Observed fire size') +
  ylab("Expected fire size")



# Visualize changes in exceedance sigma through time ----------------------
# burn areas
sigma_mu_df <- mu_df %>%
  filter(response == 2) %>%
  full_join(tbl_df(burn_covs)) %>%
  left_join(dplyr::distinct(tbl_df(ecoregions),
                            NA_L3NAME, NA_L1CODE)) %>%
  mutate(facet_factor = paste0(NA_L1CODE, NA_L3NAME))

sigma_mu_df %>%
  ggplot(aes(ym, exp(mode))) +
  geom_ribbon(aes(ymin = exp(lo), ymax = exp(hi)),
              alpha = .5, fill = 'dodgerblue',
              color = NA) +
  geom_line(size = .1) +
  facet_wrap(~ NA_L3NAME) +
  xlab('Date') +
  ylab("Standard deviation in burn area log exceedance") +
  scale_y_log10()




# Visualize covariate effects on exceedance -------------------------------
plot_size_vs_var <- function(burn_mu_df, var) {
  burn_mu_df %>%
    right_join(as_tibble(train_burns)) %>%
    distinct(ym, NA_L3NAME, .keep_all = TRUE) %>%
    ggplot(aes_string(var, "exp(mode)",
                      color = "NA_L1NAME")) +
    theme_minimal() +
    geom_segment(aes_string(xend = var,
                            y = "exp(lo)",
                            yend = "exp(hi)"),
                 alpha = .1) +
    geom_point(shape = 1, size = .1) +
    facet_wrap(~ paste(NA_L1CODE, NA_L2NAME, sep = ':')) +
    xlab(var) +
    ylab("Expected burn area exceedance over 1000 acres") +
    scale_y_log10() +
    theme(strip.text = element_text(size=6),
          legend.position = 'none') +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 1.5))) +
    scale_color_gdocs('Level 1 ecoregion') +
    scale_x_log10()
}

pw <- 10
ph <- 6

plot_size_vs_var(burn_mu_df, var = 'pr') +
  xlab('Total precipitation (same month)')
ggsave(filename = 'fig/fire-size-precip.pdf', width = pw, height = ph)

plot_size_vs_var(burn_mu_df, var = 'prev_12mo_precip') +
  xlab('Previous 12 months precipitation')
ggsave(filename = 'fig/fire-size-12mo-precip.pdf', width = pw, height = ph)

plot_size_vs_var(burn_mu_df, var = 'tmmx') +
  xlab('Mean daily maximum air temperature')
ggsave(filename = 'fig/fire-size-tmmx.pdf', width = pw, height = ph)

plot_size_vs_var(burn_mu_df, var = 'vs') +
  xlab('Mean daily wind speed')
ggsave(filename = 'fig/fire-size-wind-speed.pdf', width = pw, height = ph)

plot_size_vs_var(burn_mu_df, var = 'pet') +
  xlab('Mean potential evapotranspiration')
ggsave(filename = 'fig/fire-size-pet.pdf', width = pw, height = ph)







# Visualize effects on expected number of counts --------------------------
train_counts$row_c <- 1:nrow(train_counts)
c_df <- mu_df %>%
  filter(response == 3) %>%
  full_join(as_tibble(burn_covs))

c_df <- c_df %>%
  left_join(dplyr::distinct(tbl_df(ecoregions),
                            NA_L3NAME, NA_L1CODE)) %>%
  mutate(facet_factor = paste0(NA_L1CODE, NA_L3NAME))

plot_c_ts <- function(df) {
  df  %>%
    ggplot(aes(x=ym, y = exp(mode))) +
    theme_minimal() +
    geom_ribbon(aes(ymin = exp(lo), ymax = exp(hi),
                    fill = NA_L2NAME),
                alpha = .8,
                col = NA) +
    geom_line(size = .2) +
    facet_wrap(~facet_factor) +
    xlab('Time') +
    ylab('Expected # of burns over 1000 acres') +
    scale_y_log10() +
    theme(legend.position = 'none') +
    geom_vline(xintercept = cutoff_year, linetype = 'dashed', col = 'grey') +
    scale_fill_gdocs()
}

c_df %>%
  plot_c_ts

c_df %>%
  filter(NA_L1NAME == 'EASTERN TEMPERATE FORESTS') %>%
  plot_c_ts

c_df %>%
  filter(NA_L1NAME == 'GREAT PLAINS') %>%
  plot_c_ts

c_df %>%
  filter(NA_L1NAME == 'NORTH AMERICAN DESERTS') %>%
  plot_c_ts

c_df %>%
  filter(NA_L1NAME == 'NORTHWESTERN FORESTED MOUNTAINS') %>%
  plot_c_ts

# plot some on natural scale
c_df %>%
  filter(NA_L3NAME == 'Southern Rockies') %>%
  plot_c_ts +
  scale_y_continuous()

c_df %>%
  filter(NA_L3NAME == 'Canadian Rockies') %>%
  plot_c_ts +
  scale_y_continuous()

c_df %>%
  filter(NA_L3NAME == 'Idaho Batholith') %>%
  plot_c_ts +
  scale_y_continuous()

c_df %>%
  filter(NA_L3NAME == 'Piedmont') %>%
  plot_c_ts +
  scale_y_continuous()

c_df %>%
  filter(NA_L3NAME == 'Cross Timbers') %>%
  plot_c_ts +
  scale_y_continuous()



plot_c_vs_var <- function(df, var) {
  df %>%
    distinct(ym, NA_L3NAME, .keep_all = TRUE) %>%
    filter(year < cutoff_year) %>%
    ggplot(aes_string(var, "exp(mode - log(area))",
                      color = "NA_L1NAME")) +
    scale_x_log10() +
    theme_minimal() +
    geom_segment(aes_string(xend = var,
                            y = "exp(lo - log(area))",
                            yend = "exp(hi - log(area))"),
                 alpha = .02) +
    geom_point(shape = 1, size = .1, alpha = .3) +
    facet_wrap(~ paste(NA_L1CODE, NA_L2NAME, sep = ':')) +
    xlab(var) +
    ylab("Expected # fires > 1000 acres per sq. meter") +
    scale_y_log10() +
    theme(strip.text = element_text(size=6),
          legend.position = 'none') +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 1.5))) +
    scale_color_gdocs()
}

plot_c_vs_var(c_df, 'pr') +
  xlab('Monthly precipitation')
ggsave(filename = 'fig/fire-num-precip.pdf', width = pw, height = ph)

plot_c_vs_var(c_df, var = 'prev_12mo_precip') +
  xlab('Previous 12 months precipitation')
ggsave(filename = 'fig/fire-num-12mo-precip.pdf', width = pw, height = ph)

plot_c_vs_var(c_df, 'tmmx') +
  xlab('Mean daily maximum air temperature')
ggsave(filename = 'fig/fire-num-tmmx.pdf', width = pw, height = ph)

plot_c_vs_var(c_df, 'vs') +
  xlab('Mean daily wind speed')
ggsave(filename = 'fig/fire-num-wind-speed.pdf', width = pw, height = ph)

plot_c_vs_var(c_df, 'pet') +
  xlab('Mean daily potential evapotranspiration')
ggsave(filename = 'fig/fire-num-pet.pdf', width = pw, height = ph)


# Mapping covariate effects -----------------------------------------------

# to map the effect of variable x for an l3 ecoregion we need:
# main effect of x
# + l1_adj
# + l2_adj
# + l3_adj
beta_df <- beta_df %>%
  mutate(variable = colnames(X)[col])

er_names <- distinct(ecoregion_df, NA_L3NAME, NA_L2NAME, NA_L1NAME) %>%
  tbl_df

effect_combos <- expand.grid(NA_L3NAME = er_names$NA_L3NAME,
                             vars = c('cpr', 'ctmx',
                                      'cvs', 'cpr12',
                                      'cyear')) %>%
  tbl_df %>%
  left_join(er_names)
effect_combos$beta <- NA
effect_combos$beta_c <- NA

for (i in 1:nrow(effect_combos)) {
  var_idx <- grepl(effect_combos$vars[i], beta_df$variable)
  main_eff_idx <- grepl(paste0('^', effect_combos$vars[i], '$'),
                        beta_df$variable)
  index <-  main_eff_idx +
    var_idx * grepl(effect_combos$NA_L1NAME[i], beta_df$variable) +
    var_idx * grepl(effect_combos$NA_L2NAME[i], beta_df$variable) +
    var_idx * grepl(effect_combos$NA_L3NAME[i], beta_df$variable)

  effect_combos$beta[i] <- sum(beta_df$mode[index > 0 & beta_df$dim == 1])
  effect_combos$beta_c[i] <- sum(beta_df$mode[index > 0 & beta_df$dim == 3])
}

# simplify ecoregion data frame for easy plotting
simple_ecoregions <- ecoregions %>%
  as('Spatial') %>%
  ms_simplify(keep = 0.01) %>%
  as('sf')

ggplot() +
  geom_sf(data = simple_ecoregions,
          aes(fill = Shape_Area))

effect_sf <- effect_combos %>%
  gather(which_beta, posterior_mode, -vars, -starts_with('NA')) %>%
  full_join(simple_ecoregions) %>%
  mutate(response = ifelse(which_beta == 'beta',
                           'Burn area exceedance > 1000 acres', 'Number of fires > 1000 acres'),
         variable = ifelse(vars == 'cpr',
                           'Precipitation (same month)',
                           ifelse(vars == 'cpr12',
                                  'Precipitation (previous 12 months)',
                                  ifelse(vars == 'cvs',
                                         'Wind speed',
                                         ifelse(vars == 'ctmx',
                                                'Air temperature',
                                                'Time trend')))))

lowcolor <- 'royalblue4'
hicolor <- 'red3'

effect_map <- function(df) {
  ggplot(df, aes(fill = posterior_mode)) +
    geom_sf(size = .1, color = scales::alpha(1, .5)) +
    facet_wrap( ~ variable, nrow = 1, strip.position = 'bottom') +
    scale_fill_gradient2(low = lowcolor, high = hicolor, "") +
    theme_minimal()+
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_line(color = NA),
          panel.spacing = unit(0, "lines"))
}

p1 <- effect_sf %>%
  filter(which_beta == 'beta_c', vars != 'cyear') %>%
  effect_map() +
  ggtitle('A. Estimated effects: number of fires > 1000 acres')

p2 <- effect_sf %>%
  filter(which_beta == 'beta', vars != 'cyear') %>%
  effect_map() +
  ggtitle('B. Estimated effects: burn area exceedance > 1000 acres')

plot_grid(p1, p2, nrow = 2)
ggsave(filename = 'fig/climatic-effects.pdf', width = 16, height = 6)

# plot climate independent linear time trends
cplot <- effect_sf %>%
  filter(vars == 'cyear', which_beta == 'beta_c') %>%
  ggplot(aes(fill = posterior_mode)) +
  geom_sf(size = .1, color = scales::alpha(1, .5)) +
  scale_fill_gradient2(low = lowcolor, high = hicolor, "") +
  theme_minimal()+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_line(color = NA),
        panel.spacing = unit(0, "lines"))  +
  ggtitle('A. Estimated annual time trend: number of fires > 1000 acres')

splot <- effect_sf %>%
  filter(vars == 'cyear', which_beta == 'beta') %>%
  ggplot(aes(fill = posterior_mode)) +
  geom_sf(size = .1, color = scales::alpha(1, .5)) +
  scale_fill_gradient2(low = lowcolor, high = hicolor, "") +
  theme_minimal()+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_line(color = NA),
        panel.spacing = unit(0, "lines"))  +
  ggtitle('B. Estimated annual time trend: burn area exceedance > 1000 acres')

plot_grid(cplot, splot, nrow = 2)
ggsave('fig/time-trends.pdf', width = 7, height = 10)

# Compare predicted to expected counts for training data
c_df %>%
  dplyr::select(ym, NA_L3NAME, mode, lo, hi) %>%
  right_join(train_counts) %>%
  ggplot(aes(x = n_fire, y = exp(mode))) +
  geom_point(alpha = .1) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'blue') +
  geom_segment(aes(xend = n_fire, y = exp(lo), yend = exp(hi)), alpha = .1)


## Posterior predictive checks for extremes
n_draw <- dim(post$mu)[1]
n_preds <- dim(post$mu)[3]
block_maxima <- matrix(nrow = n_draw, ncol = n_preds)
area_sums <- matrix(nrow = n_draw, ncol = n_preds)
n_events <- rpois(n_draw * n_preds,
                  lambda = exp(c(post$mu[, 3, ]))) %>%
  matrix(nrow = n_draw, ncol = n_preds)

pb <- txtProgressBar(max = n_draw, style = 3)
for (i in 1:n_draw) {
  for (j in 1:n_preds) {
    if (n_events[i, j] > 0) {
      burn_areas <- 1e3 + exp(rnorm(n = n_events[i, j],
                                    mean = post$mu[i, 1, j],
                                    sd = exp(post$mu[i, 2, j])))
      block_maxima[i, j] <- max(burn_areas)
      area_sums[i, j] <- sum(burn_areas)
    }
  }
  setTxtProgressBar(pb, i)
}

ppred_df <- reshape2::melt(block_maxima,
                           varnames = c('iter', 'idx'),
                           value.name = 'max_burn_area') %>%
  bind_cols(reshape2::melt(area_sums,
            varnames = c('iter', 'idx'),
            value.name = 'total_burn_area')) %>%
  as_tibble %>%
  dplyr::select(-iter1, -idx1)

write_rds(ppred_df, 'ppred_df.rds')
