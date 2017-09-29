library(raster)
library(rstan)
library(sf)
library(rmapshaper)
library(tidyverse)

source('R/02-explore.R')


m_fit <- read_rds('m_fit.rds')


# Evaluate convergence ----------------------------------------------------
traceplot(m_fit,
          inc_warmup = TRUE)

traceplot(m_fit, pars = c('tau', 'sigma_size', 'sigma_eps', 'sigma_mu',
                          'tau_c', 'c_sq'))

plot(m_fit, pars = 'beta') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  coord_flip()
plot(m_fit, pars = 'beta_c') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  coord_flip()


# Extract posterior draws and visualize results ---------------------------
post <- rstan::extract(m_fit)
str(post)
rm(m_fit)
gc()



# Coefficients that seem to be far from zero ------------------------------
beta_df <- post$beta %>%
  reshape2::melt(varnames = c('iter', 'col')) %>%
  tbl_df %>%
  group_by(col) %>%
  summarize(med = median(value),
            lo = quantile(value, 0.05),
            hi = quantile(value, .95)) %>%
  ungroup %>%
  mutate(variable = colnames(X)[col],
         over_zero = lo < 0 & hi > 0)

beta_df %>%
  filter(!over_zero) %>%
  ggplot(aes(med, variable)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_segment(aes(x = lo, xend = hi, yend = variable))


beta_c_df <- post$beta_c %>%
  reshape2::melt(varnames = c('iter', 'col')) %>%
  tbl_df %>%
  group_by(col) %>%
  summarize(med = median(value),
            lo = quantile(value, 0.05),
            hi = quantile(value, .95)) %>%
  ungroup %>%
  mutate(variable = colnames(Xc)[col],
         over_zero = lo < 0 & hi > 0)

beta_c_df %>%
  filter(!over_zero) %>%
  ggplot(aes(med, variable)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_segment(aes(x = lo, xend = hi, yend = variable)) +
  theme_bw()


# explore the possibility of allowing correlation between beta and beta_c
beta_df %>%
  mutate(Parameter = 'beta') %>%
  full_join(mutate(beta_c_df, Parameter = 'beta_c')) %>%
  dplyr::select(col, Parameter, med) %>%
  spread(Parameter, med) %>%
  ggplot(aes(beta_c, beta)) +
  geom_point()
# doesn't seem to be evidence for correpsondence

# Visualize some predictions ----------------------------------------------

# burn areas
mu_df <- post$mu %>%
  reshape2::melt(varnames = c('iter', 'row')) %>%
  tbl_df %>%
  group_by(row) %>%
  summarize(med = median(value),
            lo = quantile(value, 0.025),
            hi = quantile(value, .975)) %>%
  ungroup %>%
  bind_cols(tbl_df(burn_covs))

mu_df <- mu_df %>%
  left_join(dplyr::distinct(tbl_df(ecoregions),
                          NA_L3NAME, NA_L1CODE)) %>%
  mutate(facet_factor = paste0(NA_L1CODE, NA_L3NAME))

# changes over time
mu_df %>%
  ggplot(aes(ym, exp(med))) +
  geom_ribbon(aes(ymin = exp(lo), ymax = exp(hi)),
              alpha = .5, fill = 'dodgerblue',
              color = NA) +
  geom_line(size = .1) +
  facet_wrap(~ NA_L3NAME) +
  xlab('Date') +
  ylab("Expected fire size") +
  scale_y_log10()

# true vs. predicted means
mu_df %>%
  right_join(tbl_df(train_burns)) %>%
  group_by(NA_L3NAME, ym) %>%
  ggplot(aes(log(R_ACRES - 1e3), med)) +
  geom_point(size = .2, alpha = .2) +
  geom_segment(aes(xend = log(R_ACRES - 1e3),
                   y = lo, yend = hi), alpha = .2) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
  facet_wrap(~ NA_L3NAME) +
  xlab('Observed fire size') +
  ylab("Expected fire size")

# is there an association with log(lambda) that we can exploit?
plot_size_vs_var <- function(mu_df, var) {
  mu_df %>%
    distinct(ym, NA_L3NAME, .keep_all = TRUE) %>%
    ggplot(aes_string(var, "exp(med)",
                      color = "NA_L1NAME")) +
    geom_segment(aes_string(xend = var,
                            y = "exp(lo)",
                            yend = "exp(hi)"),
                 alpha = .5) +
    geom_point(shape = 1, alpha = .5, size = .1) +
    facet_wrap(~ facet_factor) +
    xlab(var) +
    ylab("Expected burn area exceedance over 1000 acres") +
    scale_y_log10() +
    theme(strip.text = element_text(size=9)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 1.5))) +
    scale_color_discrete('Level 1 ecoregion')
}

# plot_size_vs_var(mu_df, var = 'pr') +
#   xlab('Total precipitation')
# ggsave(filename = 'fig/fire-size-precip.pdf', width = 20, height = 15)
#
# plot_size_vs_var(mu_df, var = 'prev_12mo_precip') +
#   xlab('Previous 12 months precipitation')
# ggsave(filename = 'fig/fire-size-12mo-precip.pdf', width = 20, height = 15)
#
# plot_size_vs_var(mu_df, var = 'tmmx') +
#   xlab('Mean daily maximum air temperature')
# ggsave(filename = 'fig/fire-size-tmmx.pdf', width = 20, height = 15)
#
# plot_size_vs_var(mu_df, var = 'vs') +
#   xlab('Mean daily wind speed')
# ggsave(filename = 'fig/fire-size-wind-speed.pdf', width = 20, height = 15)
#
# plot_size_vs_var(mu_df, var = 'pet') +
#   xlab('Mean potential evapotranspiration')
# ggsave(filename = 'fig/fire-size-pet.pdf', width = 20, height = 15)







# Visualize effects on expected number of counts --------------------------
train_counts$row_c <- 1:nrow(train_counts)
c_df <- post$mu_counts %>%
  reshape2::melt(varnames = c('iter', 'row_c')) %>%
  tbl_df %>%
  group_by(row_c) %>%
  summarize(med_c = median(value),
            lo_c = quantile(value, 0.025),
            hi_c = quantile(value, .975)) %>%
  ungroup %>%
  bind_cols(count_covs)

c_df <- c_df %>%
  left_join(dplyr::distinct(tbl_df(ecoregions),
                          NA_L3NAME, NA_L1CODE)) %>%
  mutate(facet_factor = paste0(NA_L1CODE, NA_L3NAME))

plot_c_ts <- function(df) {
  df  %>%
    ggplot(aes(x=ym, y = exp(med_c))) +
    geom_ribbon(aes(ymin = exp(lo_c), ymax = exp(hi_c),
                    fill = NA_L1NAME),
                alpha = .8,
                col = NA) +
    geom_line(size = .2) +
    facet_wrap(~facet_factor) +
    xlab('Time') +
    ylab('Expected # of burns over 1000 acres') +
    scale_y_log10() +
    theme(legend.position = 'bottom')
}

# c_df %>%
#   plot_c_ts
#
# c_df %>%
#   filter(NA_L1NAME == 'EASTERN TEMPERATE FORESTS') %>%
#   plot_c_ts
#
# c_df %>%
#   filter(NA_L1NAME == 'GREAT PLAINS') %>%
#   plot_c_ts
#
# c_df %>%
#   filter(NA_L3NAME == 'Southern Rockies') %>%
#   plot_c_ts
#
# plot_c_vs_var <- function(df, var) {
#   df %>%
#     distinct(ym, NA_L3NAME, .keep_all = TRUE) %>%
#     ggplot(aes_string(var, "exp(med_c - log(area))",
#                       color = "NA_L1NAME")) +
#     geom_segment(aes_string(xend = var,
#                             y = "exp(lo_c - log(area))",
#                             yend = "exp(hi_c - log(area))"),
#                  alpha = .1) +
#     geom_point(shape = 1, alpha = .8, size = .1) +
#     facet_wrap(~ facet_factor) +
#     xlab(var) +
#     ylab("Expected # fires > 1000 acres per sq. meter") +
#     scale_y_log10() +
#     theme(strip.text = element_text(size=9)) +
#     guides(colour = guide_legend(override.aes = list(alpha = 1, size = 1.5)))
# }
#
# plot_c_vs_var(c_df, 'pr') +
#   xlab('Monthly precipitation')
# ggsave(filename = 'fig/fire-num-precip.pdf', width = 20, height = 15)
#
# plot_c_vs_var(c_df, var = 'prev_12mo_precip') +
#   xlab('Previous 12 months precipitation')
# ggsave(filename = 'fig/fire-num-12mo-precip.pdf', width = 20, height = 15)
#
# plot_c_vs_var(c_df, 'tmmx') +
#   xlab('Mean daily maximum air temperature')
# ggsave(filename = 'fig/fire-num-tmmx.pdf', width = 20, height = 15)
#
# plot_c_vs_var(c_df, 'vs') +
#   xlab('Mean daily wind speed')
# ggsave(filename = 'fig/fire-num-wind-speed.pdf', width = 20, height = 15)
#
# plot_c_vs_var(c_df, 'pet') +
#   xlab('Mean daily potential evapotranspiration')
# ggsave(filename = 'fig/fire-num-pet.pdf', width = 20, height = 15)


# Mapping covariate effects -----------------------------------------------

# to map the effect of variable x for an l3 ecoregion we need:
# main effect of x
# + l1_adj
# + l2_adj
# + l3_adj
colnames(post$beta) <- colnames(X)
colnames(post$beta_c) <- colnames(Xc)

beta_meds <- apply(post$beta, 2, median)
beta_c_meds <- apply(post$beta_c, 2, median)

stopifnot(all(identical(colnames(X), colnames(Xc))))

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
  var_idx <- grepl(effect_combos$vars[i], names(beta_meds))
  main_eff_idx <- grepl(paste0('^', effect_combos$vars[i], '$'),
                        names(beta_meds))
  index <-  main_eff_idx +
    var_idx * grepl(effect_combos$NA_L1NAME[i], names(beta_meds)) +
    var_idx * grepl(effect_combos$NA_L2NAME[i], names(beta_meds)) +
    var_idx * grepl(effect_combos$NA_L3NAME[i], names(beta_meds))

  effect_combos$beta[i] <- sum(beta_meds[index > 0])
  effect_combos$beta_c[i] <- sum(beta_c_meds[index > 0])
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
  gather(which_beta, posterior_med, -vars, -starts_with('NA')) %>%
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
  ggplot(df, aes(fill = posterior_med)) +
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
  ggplot(aes(fill = posterior_med)) +
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
  ggplot(aes(fill = posterior_med)) +
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
  dplyr::select(-row_c) %>%
  right_join(train_counts) %>%
  ggplot(aes(x = n_fire, y = exp(med_c))) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
  geom_segment(aes(xend = n_fire, y = exp(lo_c), yend = exp(hi_c))) +
  facet_wrap(~ NA_L3NAME)



# Posterior predictions  --------------------------------------------------

