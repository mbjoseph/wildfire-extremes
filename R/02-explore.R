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

