library(tidyverse)
library(lubridate)
library(rstan)
library(splines)
library(spdep)

source('R/merge-data.R')

# generate spatial neighbors
if (!file.exists('nb.rds')) {
  nb <- poly2nb(as(ecoregions, 'Spatial'))
  write_rds(nb, 'nb.rds')
} else {
  nb <- read_rds('nb.rds')
}

nb_agg <- aggregate(nb, ecoregions$NA_L3NAME)

nb_mat <- nb2mat(nb_agg, style = 'B')

# generate neighborhood data for car prior
listw <- nb2listw(nb_agg, style = 'B', zero.policy = TRUE)
listw$style
B <- as(listw, 'symmetricMatrix')
# B is suitable for building N, N_edges, node1, and node2
# following http://mc-stan.org/users/documentation/case-studies/icar_stan.html





ecoregion_df <- as(ecoregions, "Spatial") %>%
  data.frame

# get areas for each L3 ecoregion
area_df <- ecoregion_df %>%
  as.data.frame %>%
  tbl_df %>%
  group_by(NA_L3NAME) %>%
  summarize(area = sum(Shape_Area))

count_df <- count_df %>%
  left_join(area_df) %>%
  arrange(ym, NA_L3NAME)

er_df <- dplyr::distinct(data.frame(ecoregions),
                       NA_L3NAME, NA_L2NAME, NA_L1NAME,
                       NA_L2CODE, NA_L1CODE) %>%
  as_tibble %>%
  filter(NA_L2NAME != 'UPPER GILA MOUNTAINS (?)')

st_covs <- ecoregion_summaries %>%
  left_join(er_df) %>%
  filter(!NA_L2NAME == "UPPER GILA MOUNTAINS (?)",
         year > 1983,
         ym <= max(mtbs$ym)) %>%
  mutate(log_housing_density = log(housing_density),
         pr = ifelse(pr < 0 , 0, pr)) %>%
  left_join(area_df) %>%
  droplevels %>%
  mutate(er_ym = paste(NA_L3NAME, ym, sep = "_")) %>%
  arrange(ym, NA_L3NAME)


assert_that(length(setdiff(st_covs$NA_L3NAME, count_df$NA_L3NAME)) == 0)
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

holdout_counts <- count_df %>%
  filter(year >= cutoff_year)

holdout_c_idx <- match(holdout_counts$er_ym, st_covs$er_ym)

holdout_burns <-  mtbs %>%
  filter(FIRE_YEAR >= cutoff_year) %>%
  left_join(st_covs)

holdout_b_idx <- match(holdout_burns$er_ym, holdout_counts$er_ym)


# this data frame has no duplicate ecoregion X timestep combos
N <- length(unique(st_covs$NA_L3NAME))
T <- length(unique(st_covs$ym))

assert_that(identical(nrow(st_covs), N * T))

# Create b-splines for climate vars
vars <- c('log_housing_density', 'vs',
          'pr', 'prev_12mo_precip', 'tmmx',
          'rmin')

df_each <- 5
X_bs <- list()
X_bs_df <- list()
for (i in seq_along(vars)) {
  X_bs[[i]] <- bs(st_covs[[vars[i]]], df = df_each, intercept = TRUE)

  X_bs_df[[i]] <- X_bs[[i]] %>%
    as_tibble
  names(X_bs_df[[i]]) <- paste('bs', vars[[i]], 1:df_each, sep = '_')
}
X_bs_df <- bind_cols(X_bs_df)
assert_that(!any(is.na(X_bs_df)))

# Create design matrices --------------------------------------------------
l3_terms <- paste0('NA_L3NAME * ', names(X_bs_df)) %>%
  paste(collapse = ' + ')
l2_terms <- paste0('NA_L2NAME * ', names(X_bs_df)) %>%
  paste(collapse = ' + ')
l1_terms <- paste0('NA_L1NAME * ', names(X_bs_df)) %>%
  paste(collapse = ' + ')

st_covs <- st_covs %>%
  bind_cols(lapply(X_bs_df, c)) %>%
  as_tibble

X <- model.matrix(as.formula(paste('~ 0 + ',
                                   l1_terms,
                                   l2_terms,
                                   l3_terms,
                                   sep = ' + ')),
      data = st_covs)

sparse_X <- extract_sparse_parts(X)
colnamesX <- colnames(X)


# design matrix for training counts
# is a subset of the rows of X, based on which rows show up in train_counts
eps_idx_train <- match(train_counts$er_ym, st_covs$er_ym)
X_tc <- X[eps_idx_train, ]
assert_that(identical(nrow(X_tc), nrow(train_counts)))
sparse_X_tc <- extract_sparse_parts(X_tc)

eps_idx_future <- setdiff(1:nrow(st_covs), eps_idx_train)
assert_that(all(diff(eps_idx_train) == 1))
assert_that(all(diff(eps_idx_future) == 1))
assert_that(eps_idx_train[length(eps_idx_train)] + 1 == eps_idx_future[1])

# design matrix for training burn areas
# is a subset of X, based on which unique rows are in train_burns
train_burn_covs <- train_burns %>%
  distinct(er_ym, .keep_all = TRUE)

# train_burn_covs has no duplicate er_ym's: should be fewer rows than train_burns
assert_that(nrow(train_burn_covs) < nrow(train_burns))
tb_idx <- match(train_burn_covs$er_ym, st_covs$er_ym)
X_tb <- X[tb_idx, ]
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

#rm(X)
rm(X_tc)
rm(X_tb)
gc()

assert_that(identical(unique(area_df$NA_L3NAME),
                      unique(st_covs$NA_L3NAME)))
assert_that(identical(levels(factor(area_df$NA_L3NAME)),
                      levels(factor(st_covs$NA_L3NAME))))



# Bundle up data into a list too pass to Stan -----------------------------
min_size <- 1e3
stan_d <- list(
  N = N,
  T = T,
  p = length(colnamesX),
  
  n_count = nrow(train_counts),
  counts = train_counts$n_fire,
  
  log_area = log(area_df$area * 1e-11) / 2,
  er_idx_train = as.numeric(factor(train_counts$NA_L3NAME,
                                   levels = levels(factor(area_df$NA_L3NAME)))),
  er_idx_full = as.numeric(factor(st_covs$NA_L3NAME)),
  
  n_fire = nrow(train_burns),
  sizes = train_burns$Acres - min_size,
  burn_idx = burn_idx,
  
  n_w = length(sparse_X$w),
  w = sparse_X$w,
  v = sparse_X$v,
  u = sparse_X$u,
  
  # sparse design matrix for training counts
  n_w_tc = length(sparse_X_tc$w),
  w_tc = sparse_X_tc$w,
  v_tc = sparse_X_tc$v,
  n_u_tc = length(sparse_X_tc$u),
  u_tc = sparse_X_tc$u,
  
  # sparse design matrix for training burns
  n_w_tb = length(sparse_X_tb$w),
  w_tb = sparse_X_tb$w,
  v_tb = sparse_X_tb$v,
  n_u_tb = length(sparse_X_tb$u),
  u_tb = sparse_X_tb$u,
  
  burn_eps_idx = burn_eps_idx,
  
  M = 1,
  slab_df = 5,
  slab_scale = 2,
  
  eps_idx_train = eps_idx_train,
  eps_idx_future = eps_idx_future,
  size_threshold = 0,
  
  n_holdout_c = length(holdout_c_idx),
  holdout_c_idx = holdout_c_idx,
  holdout_c = holdout_counts$n_fire,
  
  n_holdout_b = length(holdout_b_idx),
  holdout_b_idx = holdout_b_idx,
  holdout_b = holdout_burns$Acres - min_size,
  
  min_size = min_size,
  
  n_edges = length(B@i),
  node1 = B@i + 1, # add one to offset zero-based index
  node2 = B@j + 1,
  tb_idx = tb_idx, 
  cutoff_year = cutoff_year
)

# assert that there are no missing values in stan_d
assert_that(!any(lapply(stan_d, function(x) any(is.na(x))) %>% unlist))

zi_d <- stan_d
zi_d$M <- 2
write_rds(zi_d, 'data/processed/zi_d.rds')

write_rds(stan_d, path = 'data/processed/stan_d.rds')
print('stan_d.rds written!')


write_rds(st_covs, 'data/processed/st_covs.rds')
write_rds(cutoff_year, 'data/processed/cutoff_year.rds')
write_rds(train_counts, 'data/processed/train_counts.rds')
write_rds(holdout_counts, 'data/processed/holdout_counts.rds')
write_rds(train_burns, 'data/processed/train_burns.rds')
write_rds(holdout_burns, 'data/processed/holdout_burns.rds')
write_rds(colnamesX, 'data/processed/colnamesX.rds')
write_rds(X, 'data/processed/X.rds')
write_rds(ecoregions, 'data/processed/ecoregions.rds')
write_rds(vars, 'data/processed/vars.rds')
write_rds(mtbs, 'data/processed/mtbs.rds')
