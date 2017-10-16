library(raster)
library(rstan)
library(sf)
library(rmapshaper)
library(tidyverse)
library(magick)
library(snowfall)

source('R/02-explore.R')

most_recent_fit <- list.files(pattern = 'wfit')

m_fit <- read_rds(most_recent_fit)


# Evaluate convergence ----------------------------------------------------
traceplot(m_fit, inc_warmup = TRUE)

traceplot(m_fit, pars = c('tau', 'sigma', 'c', 'alpha'))
traceplot(m_fit, pars = c('Rho_beta', 'Rho_eps'))
plot(m_fit, pars = 'beta') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  coord_flip()


# Extract posterior draws and visualize results ---------------------------
post <- rstan::extract(m_fit)
str(post)
rm(m_fit)
gc()




# Visualize some predictions ----------------------------------------------
st_covs$row <- 1:nrow(st_covs)

mu_df <- post$mu_full %>%
  reshape2::melt(varnames = c('iter', 'response', 'row')) %>%
  tbl_df %>%
  group_by(response, row) %>%
  summarize(median = median(value),
            lo = quantile(value, 0.05),
            hi = quantile(value, 0.95)) %>%
  ungroup

## Visualize expected burn area exceedance values
burn_mu_df <- mu_df %>%
  filter(response == 1) %>%
  full_join(st_covs) %>%
  left_join(dplyr::distinct(tbl_df(ecoregions),
                            NA_L3NAME, NA_L1CODE)) %>%
  mutate(facet_factor = paste0(NA_L1CODE, NA_L3NAME))

# changes over time
plot_mu_ts <- function(df) {
  df %>%
    ggplot(aes(ym, exp(median), fill = NA_L1NAME)) +
    theme_minimal() +
    geom_ribbon(aes(ymin = exp(lo), ymax = exp(hi)),
                alpha = .5, color = NA) +
    geom_line(size = .1) +
    facet_wrap(~ NA_L3NAME) +
    xlab('Date') +
    ylab("Expected burn area log exceedance") +
    scale_y_log10() +
    geom_vline(xintercept = cutoff_year,
               linetype = 'dashed', col = 'grey') +
    scale_color_gdocs() +
    scale_fill_gdocs() +
    theme(legend.position = 'none')
}

unique(burn_mu_df$NA_L1NAME)

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
  ggplot(aes(log(R_ACRES - 1e3), median)) +
  geom_point(size = .2, alpha = .2) +
  geom_segment(aes(xend = log(R_ACRES - 1e3),
                   y = lo, yend = hi), alpha = .2) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
  facet_wrap(~ NA_L3NAME) +
  xlab('Observed fire size') +
  ylab("Expected fire size")



# Visualize changes in exceedance sigma through time ----------------------
# burn areas
mu_df %>%
  filter(response == 2) %>%
  full_join(st_covs) %>%
  left_join(dplyr::distinct(tbl_df(ecoregions),
                            NA_L3NAME, NA_L1CODE)) %>%
  mutate(facet_factor = paste0(NA_L1CODE, NA_L3NAME)) %>%
  ggplot(aes(ym, exp(median))) +
  geom_ribbon(aes(ymin = exp(lo), ymax = exp(hi)),
              alpha = .5, fill = 'dodgerblue',
              color = NA) +
  geom_line(size = .1) +
  facet_wrap(~ NA_L3NAME) +
  xlab('Date') +
  ylab("Standard deviation in burn area log exceedance") +
  scale_y_log10()




# Visualize covariate effects on exceedance -------------------------------
label_size <- 7
plot_size_vs_var <- function(burn_mu_df, var) {
  burn_mu_df %>%
    right_join(as_tibble(train_burns)) %>%
    distinct(ym, NA_L3NAME, .keep_all = TRUE) %>%
    ggplot(aes_string(var, "exp(median)",
                      color = "NA_L1NAME")) +
    theme_minimal() +
    geom_segment(aes_string(xend = var,
                            y = "exp(lo)",
                            yend = "exp(hi)"),
                 alpha = .1) +
    geom_point(shape = 1, size = .3) +
    facet_wrap(~ paste(NA_L1CODE, NA_L3NAME, sep = ':')) +
    xlab(var) +
    ylab("Expected burn area exceedance over 1000 acres") +
    scale_y_log10() +
    theme(strip.text = element_text(size=label_size),
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

plot_size_vs_var(burn_mu_df, var = 'rmin') +
  xlab('Mean minimum relative humidity')
ggsave(filename = 'fig/fire-size-rmin.pdf', width = pw, height = ph)

plot_size_vs_var(burn_mu_df, var = 'housing_density') +
  xlab('Mean housing density')
ggsave(filename = 'fig/fire-size-housing-den.pdf', width = pw, height = ph)





# Visualize effects on expected number of counts --------------------------
train_counts$row_c <- 1:nrow(train_counts)
c_df <- mu_df %>%
  filter(response == 3) %>%
  full_join(st_covs)

c_df <- c_df %>%
  left_join(dplyr::distinct(tbl_df(ecoregions),
                            NA_L3NAME, NA_L2CODE, NA_L1CODE)) %>%
  mutate(facet_factor = paste0(NA_L2CODE, NA_L3NAME))

plot_c_ts <- function(df) {
  df  %>%
    ggplot(aes(x=ym, y = exp(median))) +
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
    ggplot(aes_string(var, "exp(median - log(area))",
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
    theme(strip.text = element_text(size=label_size),
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

plot_c_vs_var(c_df, 'rmin') +
  xlab('Mean daily minimum relative humidity')
ggsave(filename = 'fig/fire-num-pet.pdf', width = pw, height = ph)

plot_c_vs_var(c_df, 'housing_density') +
  xlab('Mean housing density') +
  facet_wrap(~ NA_L3NAME, scales = 'free')
ggsave(filename = 'fig/fire-num-pet.pdf', width = pw, height = ph)


# Mapping covariate effects -----------------------------------------------

# to map the effect of variable x for an l3 ecoregion we need:
# main effect of x
# + l1_adj
# + l2_adj
# + l3_adj
# Coefficients that seem to be far from zero ------------------------------
beta_df <- post$beta %>%
  reshape2::melt(varnames = c('iter', 'dim', 'col')) %>%
  tbl_df %>%
  group_by(dim, col) %>%
  summarize(median = median(value),
            lo = quantile(value, 0.05),
            hi = quantile(value, 0.95),
            p_neg = mean(value < 0),
            p_pos = mean(value > 0)) %>%
  ungroup %>%
  mutate(variable = colnamesX[col],
         nonzero = p_neg > .9 | p_pos > .9)

beta_df %>%
  filter(nonzero) %>%
  ggplot(aes(median, variable)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_segment(aes(x = lo, xend = hi, yend = variable)) +
  facet_wrap(~ dim, scales = 'free') +
  theme(axis.text.y = element_text(size = 7))

beta_df %>%
  ggplot(aes(median, variable)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_segment(aes(x = lo, xend = hi, yend = variable)) +
  facet_wrap(~ dim, scales = 'free_x')



beta_df %>%
  select(-lo, -hi, -nonzero, -p_neg, -p_pos) %>%
  spread(dim, median) %>%
  ggplot(aes(`1`, `3`)) +
  geom_point(alpha = .5) +
  geom_vline(xintercept = 0,linetype = 'dashed') +
  geom_hline(yintercept = 0,linetype = 'dashed') +
  xlab('Effect on number of events') +
  ylab('Effect on event size')



search_string <- paste(
  c('ctmx', 'cpr', 'cpr12', 'chd', 'crmin',
    'cvs', 'cyear'), # excludes tri, since effect is constant
  collapse = '|'
)

l3_slopes <- beta_df %>%
  filter(grepl('NA_L3', variable),
         grepl(':', variable)) %>%
  select(dim, col, median, variable) %>%
  mutate(NA_L3NAME = gsub(search_string, '', variable),
         NA_L3NAME = gsub(':', '', NA_L3NAME),
         NA_L3NAME = gsub('NA_L3NAME', '', NA_L3NAME)) %>%
  group_by(dim, col) %>%
  mutate(which_coef = gsub(NA_L3NAME, '', variable),
         which_coef = gsub('NA_L3NAME:|:NA_L3NAME', '', which_coef)) %>%
  ungroup %>%
  rename(l3med = median) %>%
  select(-variable, - col)

# create dataframe will all combos of ecoregion, coef, and dim
combo_df <- expand.grid(NA_L3NAME = unique(er_df$NA_L3NAME),
                        dim = unique(l3_slopes$dim),
                        which_coef = unique(l3_slopes$which_coef)) %>%
  as_tibble %>%
  arrange(NA_L3NAME)

l3_slopes <- full_join(l3_slopes, combo_df) %>%
  arrange(NA_L3NAME, dim, which_coef) %>%
  full_join(er_df)

l2_slopes <- beta_df %>%
  filter(grepl('NA_L2', variable),
         grepl(':', variable)) %>%
  select(dim, col, median, variable)%>%
  mutate(NA_L2NAME = gsub(search_string, '', variable),
         NA_L2NAME = gsub(':', '', NA_L2NAME),
         NA_L2NAME = gsub('NA_L2NAME', '', NA_L2NAME)) %>%
  group_by(dim, col) %>%
  mutate(which_coef = gsub(NA_L2NAME, '', variable),
         which_coef = gsub('NA_L2NAME:|:NA_L2NAME', '', which_coef)) %>%
  ungroup %>%
  rename(l2med = median) %>%
  select(-variable, - col)

l2_in <- unique(er_df$NA_L2NAME) %in% l2_slopes$NA_L2NAME
assert_that(unique(er_df$NA_L2NAME)[!l2_in] ==
              unique(er_df$NA_L2NAME) %>% sort %>% `[`(1))


l1_slopes <- beta_df %>%
  filter(grepl('NA_L1', variable),
         grepl(':', variable)) %>%
  select(dim, col, median, variable) %>%
  mutate(NA_L1NAME = gsub(search_string, '', variable),
         NA_L1NAME = gsub(':', '', NA_L1NAME),
         NA_L1NAME = gsub('NA_L1NAME', '', NA_L1NAME)) %>%
  group_by(dim, col) %>%
  mutate(which_coef = gsub(NA_L1NAME, '', variable),
         which_coef = gsub('NA_L1NAME:|:NA_L1NAME', '', which_coef)) %>%
  ungroup %>%
  rename(l1med = median) %>%
  select(-variable, - col)

l1_in <- unique(er_df$NA_L1NAME) %in% l1_slopes$NA_L1NAME
assert_that(unique(er_df$NA_L1NAME)[!l1_in] ==
              unique(er_df$NA_L1NAME) %>% sort %>% `[`(1))

overall_effs <- beta_df %>%
  filter(!grepl('NA_', variable)) %>%
  select(dim, median, variable) %>%
  rename(which_coef = variable)

effect_combos <- l3_slopes %>%
  left_join(l2_slopes) %>%
  left_join(l1_slopes) %>%
  left_join(overall_effs) %>%
  filter(NA_L2NAME != 'UPPER GILA MOUNTAINS (?)') %>%
  arrange(NA_L3NAME) %>%
  group_by(dim, NA_L3NAME, which_coef) %>%
  # need to deal with R's alphabetical slope nonsense, which will have
  # NA values for the alphabetically first L1-L3 ecoregions
  summarize(l3_slope = sum(median, l1med, l2med, l3med, na.rm = TRUE)) %>%
  ungroup

# simplify ecoregion data frame for easy plotting
simple_ecoregions <- ecoregions %>%
  as('Spatial') %>%
  ms_simplify(keep = 0.005) %>%
  as('sf')

ggplot() +
  geom_sf(data = simple_ecoregions,
          aes(fill = Shape_Area))

effect_sf <- effect_combos %>%
  full_join(simple_ecoregions) %>%
  mutate(response = case_when(
    .$dim == 1 ~ 'Expected burn area exceedance > 1000 acres',
    .$dim == 2 ~ 'St. dev. burn area exceedance > 1000 acres',
    .$dim == 3 ~ 'Number of fires > 1000 acres'),
    variable = case_when(
      .$which_coef == 'cpr' ~ 'Main effect: precipitation (same month)',
      .$which_coef == 'cpr12' ~ 'Main effect: precipitation (prev. 12 months)',
      .$which_coef == 'cvs' ~ 'Main effect: wind speed',
      .$which_coef == 'ctmx' ~ 'Main effect: air temperature',
      .$which_coef == 'cyear' ~ 'Main effect: time trend',
      .$which_coef == 'chd' ~ 'Main effect: housing density',
      .$which_coef == 'crmin' ~ 'Main effect: minimum relative humidity',
      .$which_coef == 'ctmx:cvs' ~ '2 way interaction: wind speed & air temp.',
      .$which_coef == 'crmin:ctmx' ~ '2 way interaction: humidity & air temp.',
      .$which_coef == 'crmin:cvs' ~ '2 way interaction: humidity & wind speed',
      .$which_coef == 'crmin:ctmx:cvs' ~ '3 way interaction: humidity, air temp. and wind speed')
  )

lowcolor <- 'royalblue4'
hicolor <- 'red3'

effect_map <- function(df) {
  ggplot(df, aes(fill = l3_slope)) +
    geom_sf(size = .1, color = scales::alpha(1, .5)) +
    facet_wrap( ~ which_coef, nrow = 1, strip.position = 'bottom') +
    scale_fill_gradient2(low = lowcolor, high = hicolor, "") +
    theme_minimal()+
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_line(color = NA),
          panel.spacing = unit(0, "lines"))
}

p1 <- effect_sf %>%
  filter(dim == 1) %>%
  effect_map() +
  ggtitle('A. Estimated effects: expected burn area exceedance > 1000 acres')

p2 <- effect_sf %>%
  filter(dim == 2) %>%
  effect_map() +
  ggtitle('B. Estimated effects: variation in burn area exceedance > 1000 acres')

p3 <- effect_sf %>%
  filter(dim == 3) %>%
  effect_map() +
  ggtitle('C. Estimated effects: number of fires > 1000 acres')

plot_grid(p1, p2, p3, ncol = 1)
ggsave(filename = 'fig/climatic-effects.pdf', width = 16, height = 9)



# Compare predicted to expected counts for training data
c_df %>%
  dplyr::select(ym, NA_L3NAME, median, lo, hi) %>%
  right_join(count_df) %>%
  full_join(er_df) %>%
  mutate(is_train = year < cutoff_year) %>%
  filter(!is_train) %>%
  ggplot(aes(x = n_fire, y = exp(median), color = is_train)) +
  facet_wrap(~ NA_L1NAME) +
  geom_point(alpha = .5) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'blue') +
  geom_segment(aes(xend = n_fire, y = exp(lo), yend = exp(hi)), alpha = .1) +
  scale_x_log10() +
  scale_y_log10()


## Posterior predictive checks for extremes
n_draw <- dim(post$mu_full)[1]
n_preds <- dim(post$mu_full)[3]
block_maxima <- matrix(nrow = n_draw, ncol = n_preds)
area_sums <- matrix(nrow = n_draw, ncol = n_preds)
n_events <- rpois(n_draw * n_preds,
                  lambda = exp(c(post$mu_full[, 3, ]))) %>%
  matrix(nrow = n_draw, ncol = n_preds)
sum(is.na(n_events))

pb <- txtProgressBar(max = n_draw, style = 3)
for (i in 1:n_draw) {
  for (j in 1:n_preds) {
    if (!is.na(n_events[i, j]) & n_events[i, j] > 0) {
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
  bind_cols(reshape2::melt(n_events,
                           varnames = c('iter', 'idx'),
                           value.name = 'n_events')) %>%
  as_tibble %>%
  dplyr::select(-ends_with('1'), -ends_with('2')) %>%
  filter(!is.na(n_events))

write_rds(ppred_df, 'ppred_df.rds')



# Create an animated gif showing seasonality ------------------------------
anim_df <- mu_df %>%
  select(response, row, median) %>%
  filter(response != 2) %>%
  full_join(select(st_covs, row, NA_L3NAME, ym, year, month)) %>%
  spread(response, median) %>%
  rename(mu_burn = `1`,
         log_lambda = `3`)

# create a frame for each year/month
sn_df <- anim_df %>%
  left_join(er_df)

dir.create('gif')
ym_vals <- sort(unique(anim_df$ym))

save_frame <- function(ym_val) {
  out_order <- sprintf('%03d', which(ym_vals == ym_val))
  out_name <- file.path('gif', paste('frame', out_order, sep = "-"))
  out_name <- paste0(out_name, '.png')
  frame <- sn_df %>%
    filter(ym <= ym_val) %>%
    mutate(is_current = as.numeric(ym == ym_val)) %>%
    ggplot(aes(log_lambda, mu_burn,
               alpha = is_current,
               size = is_current,
               color = month)) +
    geom_point() +
    facet_wrap(~ NA_L3NAME) +
    xlab('log(lambda)') +
    ylab('mu (fire size)') +
    ggtitle(floor(ym_val)) +
    scale_alpha(range = c(0.3, 1)) +
    scale_size(range = c(1, 4)) +
    scale_color_gradient2(midpoint = 6,
                          low = 'dodgerblue',
                          mid = 'red',
                          high = 'dodgerblue',
                          limits = c(1, 12)) +
    xlim(min(sn_df$log_lambda), max(sn_df$log_lambda)) +
    ylim(min(sn_df$mu_burn), max(sn_df$mu_burn)) +
    theme_minimal() +
    theme(legend.position = 'none')

  ggsave(filename = out_name, plot = frame, width = 18, height = 12, dpi = 60)
}

sfInit(parallel = TRUE, cpus = parallel::detectCores())
sfLibrary(tidyverse)
sfExport(list = c("ym_vals", 'sn_df'))
sfSapply(ym_vals, fun = save_frame)
sfStop()

ims <- list.files('gif', full.names = TRUE) %>%
  image_read()

animation <- image_animate(ims, fps = 25)
image_write(animation, path = 'out.gif')
