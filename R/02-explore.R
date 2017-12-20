library(tidyverse)
library(lubridate)
library(rstan)
library(splines)
source("R/01-clean_data.R")

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
  mutate(crmin = c(scale(rmin)),
         cpr = c(scale(pr)),
         ctmx = c(scale(tmmx)),
         cvs = c(scale(vs)),
         cpr12 = c(scale(prev_12mo_precip)),
         ctri = c(scale(log(tri))),
         chd = c(scale(log(housing_density)))) %>%
  left_join(area_df) %>%
  droplevels %>%
  mutate(er_ym = paste(NA_L3NAME, ym, sep = "_")) %>%
  arrange(ym, NA_L3NAME)



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
vars <- c('chd', 'cvs', 'cpr', 'cpr12', 'ctmx', 'crmin')

df_each <- 4
X_bs <- list()
for (i in seq_along(vars)) {
  X_bs[[i]] <- bs(st_covs[[vars[i]]], df = df_each, intercept = TRUE) %>%
    as_tibble
  names(X_bs[[i]]) <- paste('bs', vars[[i]], 1:df_each, sep = '_')
}
X_bs <- bind_cols(X_bs)

# Create design matrices --------------------------------------------------
hd_terms <- grep('chd', names(X_bs), value = TRUE)

model_terms <- expand.grid(hd_term = hd_terms,
            other_terms = grep('chd', names(X_bs), value = TRUE, invert = TRUE)) %>%
  as.matrix %>%
  apply(1, paste, collapse = ' * ')

l3_terms <- paste0('NA_L3NAME * ', model_terms) %>%
  paste(collapse = ' + ')
l2_terms <- paste0('NA_L2NAME * ', names(X_bs)) %>%
  paste(collapse = ' + ')
l1_terms <- paste0('NA_L1NAME * ', names(X_bs)) %>%
  paste(collapse = ' + ')

st_covs <- st_covs %>%
  bind_cols(lapply(X_bs, c)) %>%
  as_tibble

X <- model.matrix(as.formula(paste('~ 0 + ctri + ', l1_terms, l2_terms, l3_terms, sep = ' + ')),
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

#rm(X)
rm(X_tc)
rm(X_tb)
gc()
