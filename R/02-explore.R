library(scales)
library(gridExtra)
library(spdep)
library(viridis)
library(tidyverse)
library(lubridate)
library(rstan)
library(cowplot)
library(corrplot)
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
nb <- as(ecoregions, "Spatial") %>%
  poly2nb
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

er_df <- dplyr::select(data.frame(ecoregions), NA_L3NAME, NA_L2NAME, NA_L1NAME) %>%
  as_tibble

st_covs <- ecoregion_summaries %>%
  mutate(dim = 1) %>%
  full_join(mutate(ecoregion_summaries, dim = 2)) %>%
  mutate(cpet = c(scale(pet)),
         cpr = c(scale(pr)),
         ctmx = c(scale(tmmx)),
         cvs = c(scale(vs)),
         dim = factor(dim),
         timestep_factor = factor(ym),
         cyear = c(scale(year))) %>%
  left_join(er_df) %>%
  filter(!NA_L2NAME == "UPPER GILA MOUNTAINS (?)")

st_covs <- st_covs[!duplicated(st_covs), ]
st_covs$id <- 1:nrow(st_covs)

# create big sparse design matrix
X <- model.matrix(~ cpr * dim * NA_L2NAME +
                    ctmx * dim * NA_L2NAME +
                    cvs * dim * NA_L2NAME +
                    cyear * dim * NA_L2NAME,
                  data = st_covs)

N <- length(unique(count_df$NA_L3NAME))
T <- length(unique(count_df$ym))

# match counts to rows in the design matrix
assert_that(nrow(count_df) == N * T)
count_idx <- rep(NA, N * T)
for (i in 1:nrow(count_df)) {
  count_idx[i] <- which(st_covs$NA_L3NAME == count_df$NA_L3NAME[i] &
          st_covs$ym == count_df$ym[i] &
            st_covs$dim == 1)
}

# match fire events to rows in design matrix
size_idx <- rep(NA, nrow(mtbs))
for (i in 1:nrow(mtbs)) {
  size_idx[i] <- which(st_covs$NA_L3NAME == mtbs$NA_L3NAME[i] &
                         st_covs$ym == mtbs$ym[i] &
                         st_covs$dim == 2)
}



# Reduce size of X to only needed cases ----------------------------------
needed_rows <- sort(unique(c(count_idx, size_idx)))
X_small <- X[needed_rows, ]
st_covs_small <- st_covs[needed_rows, ]

for (i in 1:nrow(count_df)) {
  count_idx[i] <- which(st_covs_small$NA_L3NAME == count_df$NA_L3NAME[i] &
                          st_covs_small$ym == count_df$ym[i] &
                          st_covs_small$dim == 1)
}

for (i in 1:nrow(mtbs)) {
  size_idx[i] <- which(st_covs_small$NA_L3NAME == mtbs$NA_L3NAME[i] &
                         st_covs_small$ym == mtbs$ym[i] &
                         st_covs_small$dim == 2)
}

sparse_parts <- extract_sparse_parts(X_small)

stan_d <- list(
  N = N,
  T = T,
  L = 2,
  p = ncol(X),
  nrowX = nrow(X_small),
  n_w = length(sparse_parts$w),
  w = sparse_parts$w,
  v = sparse_parts$v,
  u = sparse_parts$u,
  counts = count_df$n_fire,
  count_idx = count_idx,
  n_fire = nrow(mtbs),
  sizes = log(mtbs$R_ACRES),
  size_idx = size_idx,
  region = as.numeric(as.factor(count_df$NA_L3NAME)),
  log_area = log(count_df$area * 1e-10)
)

m_init <- stan_model('stan/st-basis.stan')
m_fit <- sampling(m_init,
                  data = stan_d,
                  pars = c('beta', 'tau', 'sigma_size', 'mu'),
                  cores = 4,
                  init_r = .1)

traceplot(m_fit)


traceplot(m_fit, pars = c('tau', 'sigma_size'))


post <- rstan::extract(m_fit)


str(post)

st_covs_small$row <- 1:nrow(st_covs_small)

mu_df <- post$mu %>%
  reshape2::melt(varnames = c('iter', 'row')) %>%
  tbl_df %>%
  group_by(row) %>%
  summarize(med = median(value),
            lo = quantile(value, 0.025),
            hi = quantile(value, .975)) %>%
  ungroup %>%
  full_join(st_covs_small)

mu_df %>%
  filter(dim == 1, NA_L3NAME == 'Southern Rockies') %>%
  ggplot(aes(ym, med)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .5, fill = 'dodgerblue') +
  geom_line() +
  facet_wrap(~ NA_L3NAME)

mu_df %>%
  filter(dim == 1) %>%
  ggplot(aes(pr, med)) +
  geom_segment(aes(xend = pr,
                   y = lo,
                   yend = hi),
               alpha = .5) +
  geom_point(shape = 1, alpha = .5, size = .1) +
  facet_wrap(~ NA_L3NAME) +
  xlab('Total precipitation') +
  ylab("Expected number of fires")

mu_df %>%
  filter(dim == 1) %>%
  ggplot(aes(pr, med)) +
  geom_segment(aes(xend = tmmx,
                   y = lo,
                   yend = hi),
               alpha = .5) +
  geom_point(shape = 1, alpha = .5, size = .1) +
  facet_wrap(~ NA_L3NAME) +
  xlab('Mean daily maximum air temperature') +
  ylab("Expected number of fires")


mu_df %>%
  filter(dim == 1) %>%
  ggplot(aes(vs, med)) +
  geom_segment(aes(xend = vs,
                   y = lo,
                   yend = hi),
               alpha = .5) +
  geom_point(shape = 1, alpha = .5, size = .1) +
  facet_wrap(~ NA_L3NAME) +
  xlab('Mean daily wind speed') +
  ylab("Expected number of fires")



mu_df %>%
  filter(dim == 1) %>%
  ggplot(aes(pet, med)) +
  geom_segment(aes(xend = pet,
                   y = lo,
                   yend = hi),
               alpha = .5) +
  geom_point(shape = 1, alpha = .5, size = .1) +
  facet_wrap(~ NA_L3NAME) +
  xlab('Mean daily potential evapotranspiration') +
  ylab("Expected number of fires")


mu_df %>%
  filter(dim == 2) %>%
  ggplot(aes(pr, exp(med))) +
  geom_segment(aes(xend = pr,
                   y = exp(lo),
                   yend = exp(hi)),
               alpha = .5) +
  geom_point(shape = 1, alpha = .5, size = .1) +
  facet_wrap(~ NA_L3NAME) +
  xlab('Total precipitation') +
  ylab("Expected fire size (log burn area)")

mu_df %>%
  filter(dim == 2) %>%
  ggplot(aes(tmmx, med)) +
  geom_segment(aes(xend = tmmx,
                   y = lo,
                   yend = hi),
               alpha = .5) +
  geom_point(shape = 1, alpha = .5, size = .1) +
  facet_wrap(~ NA_L3NAME) +
  xlab('Mean daily maximum air temperature') +
  ylab("Expected fire size (log burn area)")


mu_df %>%
  filter(dim == 2) %>%
  ggplot(aes(vs, med)) +
  geom_segment(aes(xend = vs,
                   y = lo,
                   yend = hi),
               alpha = .5) +
  geom_point(shape = 1, alpha = .5, size = .1) +
  facet_wrap(~ NA_L3NAME) +
  xlab('Mean daily wind speed') +
  ylab("Expected fire size (log burn area)")



mu_df %>%
  filter(dim == 2) %>%
  ggplot(aes(pet, med)) +
  geom_segment(aes(xend = pet,
                   y = lo,
                   yend = hi),
               alpha = .5) +
  geom_point(shape = 1, alpha = .5, size = .1) +
  facet_wrap(~ NA_L3NAME) +
  xlab('Mean daily potential evapotranspiration') +
  ylab("Expected fire size (log burn area)")



