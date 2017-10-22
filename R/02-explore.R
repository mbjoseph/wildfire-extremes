library(scales)
library(raster)
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
library(splines)
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

ecoregion_summaries %>%
  ggplot(aes(ym, rmin)) +
  geom_line() +
  facet_wrap(~ NA_L3NAME)

er_df <- dplyr::distinct(data.frame(ecoregions),
                       NA_L3NAME, NA_L2NAME, NA_L1NAME,
                       NA_L2CODE, NA_L1CODE) %>%
  as_tibble %>%
  filter(NA_L2NAME != 'UPPER GILA MOUNTAINS (?)')

st_covs <- ecoregion_summaries %>%
  left_join(er_df) %>%
  filter(!NA_L2NAME == "UPPER GILA MOUNTAINS (?)", year > 1983) %>%
  mutate(crmin = c(scale(rmin)),
         cpr = c(scale(pr)),
         ctmx = c(scale(tmmx)),
         cvs = c(scale(vs)),
         cpr12 = c(scale(prev_12mo_precip)),
         cyear = c(scale(year)),
         ctri = c(scale(log(tri))),
         chd = c(scale(log(housing_density)))) %>%
  left_join(area_df) %>%
  droplevels %>%
  mutate(er_ym = paste(NA_L3NAME, ym, sep = "_")) %>%
  arrange(ym, NA_L3NAME)



# Compute spatial means, temporal means, and residuals
st_covs <- st_covs %>%
  group_by(NA_L3NAME) %>%
  mutate(sp_mean_crmin = mean(crmin),
         sp_mean_cpr = mean(cpr),
         sp_mean_ctmx = mean(ctmx),
         sp_mean_cvs = mean(cvs),
         sp_mean_cpr12 = mean(cpr12),
         sp_mean_chd = mean(chd)) %>%
  ungroup() %>%
  group_by(ym) %>%
  mutate(t_mean_crmin = mean(crmin),
         t_mean_cpr = mean(cpr),
         t_mean_ctmx = mean(ctmx),
         t_mean_cvs = mean(cvs),
         t_mean_cpr12 = mean(cpr12),
         t_mean_chd = mean(chd)) %>%
  ungroup() %>%
  mutate(r_crmin = crmin - t_mean_crmin - sp_mean_crmin,
         r_cpr = cpr - t_mean_cpr - sp_mean_cpr,
         r_ctmx = ctmx - t_mean_ctmx - sp_mean_ctmx,
         r_cvs = cvs - t_mean_cvs - sp_mean_cvs,
         r_cpr12 = cpr12 - t_mean_cpr12 - sp_mean_cpr12,
         r_chd = chd - t_mean_chd - sp_mean_chd)


# monthly basis functions
monthly_basis <- bs(st_covs$month, df = 5, intercept = TRUE)
colnames(monthly_basis) <- paste0('mb', 1:ncol(monthly_basis))


basis_df <- monthly_basis %>%
  as.matrix %>%
  as_tibble %>%
  lapply(FUN = c) %>%
  bind_cols

st_covs <- bind_cols(st_covs, basis_df)


assert_that(!any(duplicated(st_covs)))
st_covs$id <- 1:nrow(st_covs)

# Create training sets, including years from 1984 to cutoff_year - 1
cutoff_year <- 2010

train_counts <- count_df %>%
  filter(year < cutoff_year) %>%
  left_join(st_covs)

train_burns <- mtbs %>%
  filter(FIRE_YEAR < cutoff_year) %>%
  left_join(st_covs)

# this data frame has no duplicate ecoregion X timestep combos
N <- length(unique(st_covs$NA_L3NAME))
T <- length(unique(st_covs$ym))

assert_that(identical(nrow(st_covs), N * T))




# Create design matrices --------------------------------------------------
make_X <- function(df) {
  model.matrix(~ 0 +
                 ctri +
                 sp_mean_crmin +
                 sp_mean_cpr +
                 sp_mean_ctmx +
                 sp_mean_cvs +
                 sp_mean_cpr12 +
                 sp_mean_chd +
                 t_mean_crmin +
                 t_mean_cpr +
                 t_mean_ctmx +
                 t_mean_cvs +
                 t_mean_cpr12 +
                 t_mean_chd +
                 cyear * r_chd * r_cpr * r_crmin * r_ctmx * r_cvs * r_cpr12 * NA_L3NAME +
                 cyear * r_chd * r_cpr * r_crmin * r_ctmx * r_cvs * r_cpr12 * NA_L2NAME +
                 cyear * r_chd * r_cpr * r_crmin * r_ctmx * r_cvs * r_cpr12 * NA_L1NAME +
                 mb1 * NA_L3NAME +
                 mb1 * NA_L2NAME +
                 mb1 * NA_L1NAME +
                 mb2 * NA_L3NAME +
                 mb2 * NA_L2NAME +
                 mb2 * NA_L1NAME +
                 mb3 * NA_L3NAME +
                 mb3 * NA_L2NAME +
                 mb3 * NA_L1NAME +
                 mb4 * NA_L3NAME +
                 mb4 * NA_L2NAME +
                 mb4 * NA_L1NAME +
                 mb5 * NA_L3NAME +
                 mb5 * NA_L2NAME +
                 mb5 * NA_L1NAME,
               data = df)
}

# discard all 3-way or higher interactions
X <- make_X(st_covs)
num_interacting_variables <- lengths(regmatches(colnames(X),
                                                gregexpr(":", colnames(X))))
X <- X[, num_interacting_variables <= 1]
sparse_X <- extract_sparse_parts(X)
colnamesX <- colnames(X)


# design matrix for training counts
# is a subset of the rows of X, based on which rows show up in train_counts
eps_idx_train <- match(train_counts$er_ym, st_covs$er_ym)
X_tc <- X[eps_idx_train, ]
assert_that(identical(nrow(X_tc), nrow(train_counts)))
sparse_X_tc <- extract_sparse_parts(X_tc)

eps_idx_future <- setdiff(1:nrow(st_covs), eps_idx_train)
plot(eps_idx_train,
     xlim = c(0, nrow(st_covs)),
     ylim = c(0, nrow(st_covs)),
     type = 'l')
lines(eps_idx_future, eps_idx_future, col = 2)

# design matrix for training burn areas
# is a subset of X, based on which unique rows are in train_burns
train_burn_covs <- train_burns %>%
  distinct(er_ym, .keep_all = TRUE)

# train_burn_covs has no duplicate er_ym's: should be fewer rows than train_burns
assert_that(nrow(train_burn_covs) < nrow(train_burns))

X_tb <- X[match(train_burn_covs$er_ym, st_covs$er_ym), ]
assert_that(identical(nrow(X_tb), nrow(train_burn_covs)))
sparse_X_tb <- extract_sparse_parts(X_tb)

# indices to match epsilon parameters for burn areas to those computed for counts
burn_eps_idx <- match(train_burn_covs$er_ym, train_counts$er_ym)
assert_that(train_burn_covs$er_ym[1] == train_counts$er_ym[burn_eps_idx[1]])

# indices to match each fire event to a row in the design matrix for burns
burn_idx <- match(train_burns$er_ym, train_burn_covs$er_ym)

# check to make sure the indices were correct
assert_that(max(burn_idx) <= nrow(st_covs))
assert_that(all(train_burn_covs$NA_L3NAME[burn_idx] == train_burns$NA_L3NAME))
assert_that(all(train_burn_covs$ym[burn_idx] == train_burns$ym))

rm(X)
rm(X_tc)
rm(X_tb)
gc()

weibull_scale_adj <- 1e4
