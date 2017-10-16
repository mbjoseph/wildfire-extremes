library(raster)
library(rstan)
library(sf)
library(rmapshaper)
library(tidyverse)
library(rstan)
library(loo)

source('R/02-explore.R')


w_fit <- read_rds('w_fit.rds')

size_ll <- extract_log_lik(w_fit, parameter_name = 'loglik_f')
loo_size <- loo(size_ll)
loo_size
plot(loo_size)

count_ll <- extract_log_lik(w_fit, parameter_name = 'loglik_c')
loo_count <- loo(count_ll)
loo_count
plot(loo_count)


# Evaluate convergence ----------------------------------------------------
traceplot(w_fit, inc_warmup = TRUE)

traceplot(w_fit, pars = c('tau', 'sigma', 'c', 'alpha'))
traceplot(w_fit, pars = c('Rho_beta', 'Rho_eps'))

plot(w_fit, pars = 'beta') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  coord_flip()


# Extract posterior draws and visualize results ---------------------------
post <- rstan::extract(w_fit)
str(post)
rm(w_fit)
gc()



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
         nonzero = p_neg > .8 | p_pos > .8)

beta_df %>%
  filter(nonzero) %>%
  ggplot(aes(median, variable)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_segment(aes(x = lo, xend = hi, yend = variable)) +
  facet_wrap(~ dim, scales = 'free') +
  theme(axis.text.y = element_text(size = 7))


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
    ylab("Weibull shape parameter") +
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


# Visualize changes in scale parameter through time ----------------------
# burn areas
mu_df %>%
  filter(response == 2) %>%
  full_join(st_covs) %>%
  left_join(er_df) %>%
  mutate(facet_factor = paste0(NA_L2CODE, NA_L3NAME)) %>%
  ggplot(aes(ym, exp(median), fill = NA_L1CODE)) +
  geom_ribbon(aes(ymin = exp(lo), ymax = exp(hi)),
              alpha = .5,
              color = NA) +
  geom_line(size = .1) +
  facet_wrap(~ facet_factor) +
  xlab('Date') +
  ylab("Weibull scale parameter") +
  scale_y_log10() +
  theme(legend.position = 'none') +
  scale_fill_gdocs()



# Visualize effects on expected number of counts --------------------------
train_counts$row_c <- 1:nrow(train_counts)
c_df <- mu_df %>%
  filter(response == 3) %>%
  full_join(st_covs)

c_df <- c_df %>%
  left_join(er_df) %>%
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
  geom_segment(aes(xend = n_fire, y = exp(lo), yend = exp(hi)), alpha = .5) +
  scale_x_log10() +
  scale_y_log10(limits = c(.001, 1e4))




# Mapping covariate effects -----------------------------------------------

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
    .$dim == 1 ~ 'Weibull shape parameter',
    .$dim == 2 ~ 'Weibull scale parameter',
    .$dim == 3 ~ 'Number of fires > 1000 acres'))

lowcolor <- 'royalblue4'
hicolor <- 'red3'

effect_map <- function(df) {
  ggplot(df, aes(fill = l3_slope)) +
    geom_sf(size = .1, color = scales::alpha(1, .5)) +
    facet_wrap( ~ which_coef, strip.position = 'bottom') +
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
  ggtitle('A. Estimated effects: Weibull shape parameter')

p2 <- effect_sf %>%
  filter(dim == 2) %>%
  effect_map() +
  ggtitle('B. Estimated effects: Weibull scale parameter')

p3 <- effect_sf %>%
  filter(dim == 3) %>%
  effect_map() +
  ggtitle('C. Estimated effects: number of fires > 1000 acres')




## Posterior predictive checks ---------------------------

# set breakpoints for histogram
max_breaks <- 200
mtbs_hist <- hist(mtbs$R_ACRES, breaks = max_breaks, right = FALSE)

n_draw <- dim(post$mu_full)[1]
n_preds <- dim(post$mu_full)[3]
hist_counts <- matrix(nrow = n_draw, ncol = max_breaks*1.5)
hist_mids <- matrix(nrow = n_draw, ncol = max_breaks*1.5)
block_maxima <- matrix(nrow = n_draw, ncol = n_preds)
area_sums <- matrix(nrow = n_draw, ncol = n_preds)
n_events <- rpois(n_draw * n_preds,
                  lambda = exp(c(post$mu_full[, 3, ]))) %>%
  matrix(nrow = n_draw, ncol = n_preds)
sum(is.na(n_events))

pb <- txtProgressBar(max = n_draw, style = 3)
for (i in 1:n_draw) {
  ba_vec <- list() # to hold all fire sizes for iteration i
  counter <- 1
  for (j in 1:n_preds) {
    if (!is.na(n_events[i, j]) & n_events[i, j] > 0) {
      burn_areas <- 1e3 +
        weibull_scale_adj * rweibull(n = n_events[i, j],
                                   shape = exp(post$mu[i, 1, j]),
                                   scale = exp(post$mu[i, 2, j]))
      block_maxima[i, j] <- max(burn_areas)
      area_sums[i, j] <- sum(burn_areas)
      ba_vec[[counter]] <- burn_areas
      counter <- counter + 1
    }
  }
  ba_vec <- unlist(ba_vec)
  ba_hist <- hist(ba_vec, breaks = max_breaks, plot = FALSE, right = FALSE)
  nbins <- length(ba_hist$counts)
  hist_counts[i, 1:nbins] <- ba_hist$counts
  hist_mids[i, 1:nbins] <- ba_hist$mids
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


hist_df <- reshape2::melt(hist_counts,
                          varnames = c('iter', 'bin'),
                          value.name = 'count') %>%
  bind_cols(reshape2::melt(hist_mids,
                           varnames = c('iter', 'bin'),
                           value.name = 'midpoint')) %>%
  as_tibble %>%
  select(-ends_with('1')) %>%
  arrange(iter, bin) %>%
  na.omit %>%
  mutate(approx_burn_area = count * midpoint) %>%
  group_by(iter) %>%
  mutate(cum_burn_area = cumsum(approx_burn_area)) %>%
  ungroup

write_rds(hist_df, 'ppc_hist.rds')
write_rds(mtbs_hist, 'mtbs_hist.rds')

