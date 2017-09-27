library(scales)
library(gridExtra)
library(spdep)
library(viridis)
library(rmapshaper)
library(tidyverse)
library(lubridate)
library(rstan)
library(cowplot)
library(ggthemes)
library(sf)
source("R/01-clean_data.R")


# Visualize the distribution of number of fires
count_df %>%
  ggplot(aes(ym, n_fire)) +
  geom_line() +
  facet_wrap(~NA_L3NAME) +
  scale_y_log10()

# Visualize changes in maxima
mtbs %>%
  tbl_df %>%
  group_by(ym, NA_L3NAME) %>%
  summarize(max = max(P_ACRES)) %>%
  ggplot(aes(ym, max)) +
  geom_point() +
  facet_wrap(~ NA_L3NAME)




# Create spatial neighbors ------------------------------------------------
if (!file.exists('nb.rds')) {
  nb <- as(ecoregions, "Spatial") %>%
    poly2nb
  write_rds(nb, path = 'nb.rds')
}
nb <- read_rds('nb.rds')
W <- nb %>%
  nb2mat(zero.policy = TRUE, style = 'B')
plot(raster(W), col = viridis(10))

# aggregate neighbor matrix to ecoregion level
ecoregion_df <- as(ecoregions, "Spatial") %>%
  data.frame

# get areas for each L3 ecoregion
area_df <- ecoregion_df %>%
  as.data.frame %>%
  tbl_df %>%
  group_by(NA_L3NAME) %>%
  summarize(area = sum(Shape_Area))

count_df <- count_df %>%
  left_join(area_df)


W_ag <- spdep::aggregate.nb(nb, IDs = ecoregion_df$NA_L3NAME) %>%
  nb2mat(zero.policy = TRUE, style = "B")
plot(raster(W_ag), col = viridis(10))


# visualize spatiotemporal covariates
ecoregion_summaries %>%
  gather(var, val, -NA_L3NAME, -year, -month, -ym) %>%
  ggplot(aes(ym, val)) +
  geom_line() +
  facet_wrap(~ var, scales = "free_y") +
  scale_y_log10()

ecoregion_summaries %>%
  ggplot(aes(ym, prev_12mo_precip)) +
  geom_line() +
  facet_wrap(~ NA_L3NAME)

er_df <- dplyr::select(data.frame(ecoregions),
                       NA_L3NAME, NA_L2NAME, NA_L1NAME) %>%
  as_tibble

st_covs <- ecoregion_summaries %>%
  mutate(dim = 1) %>%
  full_join(mutate(ecoregion_summaries, dim = 2)) %>%
  mutate(cpet = c(scale(pet)),
         cpr = c(scale(pr)),
         ctmx = c(scale(tmmx)),
         cvs = c(scale(vs)),
         cpr12 = c(scale(prev_12mo_precip)),
         dim = factor(dim),
         timestep_factor = factor(ym),
         cyear = c(scale(year)),
         ctri = c(scale(log(tri)))) %>%
  left_join(er_df) %>%
  filter(!NA_L2NAME == "UPPER GILA MOUNTAINS (?)")

st_covs <- st_covs[!duplicated(st_covs), ]
st_covs$id <- 1:nrow(st_covs)

# Create training sets, including years 1984 - 2004
cutoff_year <- 2010

train_counts <- count_df %>%
  filter(year < cutoff_year) %>%
  left_join(filter(st_covs, dim == 1))

train_burns <- mtbs %>%
  filter(FIRE_YEAR < cutoff_year) %>%
  left_join(filter(st_covs, dim == 2))

# this data frame has no duplicate ecoregion X timestep combos
train_burn_covs <- st_covs %>%
  filter(dim == 2, ym >= min(train_burns$ym), ym <= max(train_burns$ym)) %>%
  distinct(NA_L3NAME, ym, .keep_all = TRUE)

N <- length(unique(train_counts$NA_L3NAME))
T <- length(unique(train_counts$ym))


# Create design matrices --------------------------------------------------
make_X <- function(df) {
  model.matrix(~ 0 +
                 ctri +
                 cyear * NA_L3NAME +
                 cyear * NA_L2NAME +
                 cyear * NA_L1NAME +
                 cpr * NA_L3NAME +
                 cpr12 * NA_L3NAME +
                 ctmx * NA_L3NAME +
                 cvs * NA_L3NAME +
                 cpr * NA_L2NAME +
                 cpr12 * NA_L2NAME +
                 ctmx * NA_L2NAME +
                 cvs * NA_L2NAME +
                 cpr * NA_L1NAME +
                 cpr12 * NA_L1NAME +
                 ctmx * NA_L1NAME +
                 cvs * NA_L1NAME,
               data = df)
}

Xc <- make_X(train_counts)
sparse_Xc <- extract_sparse_parts(Xc)

# need to ensure that all ecoregions show up here
X <- make_X(train_burn_covs)
sparse_X <- extract_sparse_parts(X)

burn_idx <- rep(NA, nrow(train_burns))
for (i in 1:nrow(train_burns)) {
  burn_idx[i] <- which(train_burn_covs$NA_L3NAME == train_burns$NA_L3NAME[i] &
                         train_burn_covs$ym == train_burns$ym[i])
}

# check to make sure the indices were correct
stopifnot(max(burn_idx) <= nrow(train_burn_covs))
stopifnot(all(train_burn_covs$NA_L3NAME[burn_idx] == train_burns$NA_L3NAME))
stopifnot(all(train_burn_covs$ym[burn_idx] == train_burns$ym))

# note whether there are columns in Xc that are not in X
colnames(Xc)[!(colnames(Xc) %in% colnames(X))]

stan_d <- list(
  N = N,
  T = T,
  L = 2,
  p = ncol(X),
  p_c = ncol(Xc),

  n_w = length(sparse_X$w),
  w = sparse_X$w,
  v = sparse_X$v,
  u = sparse_X$u,
  nrowX = nrow(X),
  n_u = length(sparse_X$u),
  burn_idx = burn_idx,

  n_wc = length(sparse_Xc$w),
  wc = sparse_Xc$w,
  vc = sparse_Xc$v,
  uc = sparse_Xc$u,

  counts = train_counts$n_fire,
  n_fire = nrow(train_burns),
  sizes = log(train_burns$R_ACRES - 1e3),
  log_area = log(train_counts$area * 1e-10)
)

m_init <- stan_model('stan/st-basis.stan')
# m_fit <- sampling(m_init,
#                   data = stan_d,
#                   pars = c('beta', 'tau', 'sigma_size',
#                            'mu', 'mu_counts',
#                            'sigma_eps',
#                            'beta_c', 'tau_c',
#                            'alpha',
#                            'c_sq'),
#                   cores = 4,
#                   init_r = .01,
#                   iter = 1000)
# write_rds(m_fit, 'm_fit.rds')
m_fit <- read_rds('m_fit.rds')

traceplot(m_fit,
          inc_warmup = TRUE)

traceplot(m_fit, pars = c('tau', 'sigma_size', 'sigma_eps', 'tau_c', 'c_sq'))

plot(m_fit, pars = 'beta') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  coord_flip()
plot(m_fit, pars = 'beta_c') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  coord_flip()


post <- rstan::extract(m_fit)
str(post)
rm(m_fit)
gc()

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
  select(col, Parameter, med) %>%
  spread(Parameter, med) %>%
  ggplot(aes(beta_c, beta)) +
  geom_point()


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
  `[`(burn_idx, ) %>%
  bind_cols(tbl_df(train_burns))

mu_df <- mu_df %>%
  left_join(dplyr::select(tbl_df(ecoregions),
                          NA_L3NAME, NA_L1CODE)) %>%
  mutate(facet_factor = paste0(NA_L1CODE, NA_L3NAME))

# changes over time
mu_df %>%
  ggplot(aes(ym, exp(med))) +
  geom_ribbon(aes(ymin = exp(lo), ymax = exp(hi)),
              alpha = .5, fill = 'dodgerblue') +
  geom_line() +
  facet_wrap(~ NA_L3NAME) +
  xlab('Date') +
  ylab("Expected fire size") +
  scale_y_log10()

# true vs. predicted means
mu_df %>%
  group_by(NA_L3NAME, ym) %>%
  summarize(emp_mean = mean(log(R_ACRES - 1e3)),
            med = unique(na.omit(med)),
            lo = unique(na.omit(lo)),
            hi = unique(na.omit(hi))) %>%
  ggplot(aes(emp_mean, med)) +
  geom_point(size = .2, alpha = .2) +
  geom_segment(aes(xend = emp_mean,
                   y = lo, yend = hi), alpha = .2) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
  facet_wrap(~ NA_L3NAME) +
  xlab('Empirical mean fire size') +
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
    ylab("Expected fire size (burn area)") +
    scale_y_log10() +
    theme(strip.text = element_text(size=9)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 1.5)))
}

plot_size_vs_var(mu_df, var = 'pr') +
  xlab('Total precipitation')
ggsave(filename = 'fig/fire-size-precip.pdf', width = 20, height = 15)

plot_size_vs_var(mu_df, var = 'prev_12mo_precip') +
  xlab('Previous 12 months precipitation')
ggsave(filename = 'fig/fire-size-12mo-precip.pdf', width = 20, height = 15)

plot_size_vs_var(mu_df, var = 'tmmx') +
  xlab('Mean daily maximum air temperature')
ggsave(filename = 'fig/fire-size-tmmx.pdf', width = 20, height = 15)

plot_size_vs_var(mu_df, var = 'vs') +
  xlab('Mean daily wind speed')
ggsave(filename = 'fig/fire-size-wind-speed.pdf', width = 20, height = 15)

plot_size_vs_var(mu_df, var = 'pet') +
  xlab('Mean potential evapotranspiration')
ggsave(filename = 'fig/fire-size-pet.pdf', width = 20, height = 15)







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
  full_join(train_counts)

c_df <- c_df %>%
  left_join(dplyr::select(tbl_df(ecoregions),
                          NA_L3NAME, NA_L1CODE)) %>%
  mutate(facet_factor = paste0(NA_L1CODE, NA_L3NAME))

plot_c_ts <- function(df) {
  df  %>%
    ggplot(aes(x=ym, y = exp(med_c))) +
    geom_ribbon(aes(ymin = exp(lo_c), ymax = exp(hi_c)),
                alpha = .5, fill = 'red', col = NA) +
    geom_line() +
    facet_wrap(~facet_factor) +
    xlab('Time') +
    ylab('Expected # of burns over 1000 acres') +
    scale_y_log10()
}

c_df %>%
  filter(NA_L1NAME == 'EASTERN TEMPERATE FORESTS') %>%
  plot_c_ts

c_df %>%
  filter(NA_L1NAME == 'GREAT PLAINS') %>%
  plot_c_ts

c_df %>%
  filter(NA_L1NAME == 'NORTHWESTERN FORESTED MOUNTAINS') %>%
  plot_c_ts

plot_c_vs_var <- function(df, var) {
  df %>%
    distinct(ym, NA_L3NAME, .keep_all = TRUE) %>%
    ggplot(aes_string(var, "exp(med_c - log(area))",
                      color = "NA_L1NAME")) +
    geom_segment(aes_string(xend = var,
                            y = "exp(lo_c - log(area))",
                            yend = "exp(hi_c - log(area))"),
                 alpha = .1) +
    geom_point(shape = 1, alpha = .8, size = .1) +
    facet_wrap(~ facet_factor) +
    xlab(var) +
    ylab("Expected # fires > 1000 acres per sq. meter") +
    scale_y_log10() +
    theme(strip.text = element_text(size=9)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 1.5)))
}

plot_c_vs_var(c_df, 'pr') +
  xlab('Monthly precipitation')
ggsave(filename = 'fig/fire-num-precip.pdf', width = 20, height = 15)

plot_c_vs_var(c_df, var = 'prev_12mo_precip') +
  xlab('Previous 12 months precipitation')
ggsave(filename = 'fig/fire-num-12mo-precip.pdf', width = 20, height = 15)

plot_c_vs_var(c_df, 'tmmx') +
  xlab('Mean daily maximum air temperature')
ggsave(filename = 'fig/fire-num-tmmx.pdf', width = 20, height = 15)

plot_c_vs_var(c_df, 'vs') +
  xlab('Mean daily wind speed')
ggsave(filename = 'fig/fire-num-wind-speed.pdf', width = 20, height = 15)

plot_c_vs_var(c_df, 'pet') +
  xlab('Mean daily potential evapotranspiration')
ggsave(filename = 'fig/fire-num-pet.pdf', width = 20, height = 15)


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
                             vars = c('cpr', 'ctmx', 'cvs', 'cpr12')) %>%
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
                           'Burn area', 'Number of fires'),
         variable = ifelse(vars == 'cpr',
                           'Precipitation (same month)',
                           ifelse(vars == 'cpr12',
                                  'Precipitation (previous 12 months)',
                                  ifelse(vars == 'cvs',
                                         'Wind speed',
                                         'Air temperature'))))

lowcolor <- 'royalblue4'
hicolor <- 'red3'

p1 <- ggplot(filter(effect_sf, which_beta == 'beta_c'),
             aes(fill = posterior_med)) +
  geom_sf(size = .1, color = scales::alpha(1, .5)) +
  facet_wrap( ~ variable, nrow = 1, strip.position = 'bottom') +
  scale_fill_gradient2(low = lowcolor, high = hicolor,
                       "") +
  theme_minimal()+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_line(color = NA),
        panel.spacing = unit(0, "lines"),
        plot.margin = margin(0, 0, 0, 0, "cm")) +
  ggtitle('Estimated climate effects on the monthly expected number of fires over 1000 acres')

p2 <- ggplot(filter(effect_sf, which_beta == 'beta'),
       aes(fill = posterior_med)) +
  geom_sf(size = .1, color = scales::alpha(1, .5)) +
  facet_wrap(~ variable, nrow = 1, strip.position = 'bottom') +
  scale_fill_gradient2(low = lowcolor, high = hicolor,
                       "") +
  theme_minimal()+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_line(color = NA),
        panel.spacing = unit(0, "lines"),
        plot.margin = margin(0, 0, 0, 0, "cm")) +
  ggtitle('Estimated climate effects on the monthly expected burn area exceedance over 1000 acres')

plot_grid(p1, p2, nrow = 2)
ggsave(filename = 'fig/climatic-effects.pdf', width = 16, height = 6)


c_df %>%
  ggplot(aes(x = n_fire, y = exp(med_c))) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
  geom_segment(aes(xend = n_fire, y = exp(lo_c), yend = exp(hi_c)))

mu_df %>%
  left_join(dplyr::select(c_df, NA_L3NAME, ym, med_c, lo_c, hi_c)) %>%
  ggplot(aes(x = med_c, y = log(R_ACRES - 1e3))) +
  geom_point(shape = 1, size = .5, alpha = .3) +
  facet_wrap(~NA_L3NAME)
