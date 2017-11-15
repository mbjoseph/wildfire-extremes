source('R/02-explore.R')
library(extraDistr)
library(ggridges)

fit <- read_rds(path = list.files(pattern = 'zinbgpdfit_.*')[3])

# Extract posterior draws and visualize results ---------------------------
post <- rstan::extract(fit)
str(post)
rm(fit)
gc()

## Posterior predictive simulations ---------------------------

# set breakpoints for histogram
max_breaks <- 200
mtbs_hist <- hist(mtbs$R_ACRES, breaks = max_breaks, right = FALSE)
write_rds(mtbs_hist, 'mtbs_hist.rds')

n_draw <- dim(post$mu_full)[1]
n_preds <- dim(post$mu_full)[3]
hist_counts <- matrix(nrow = n_draw, ncol = max_breaks*1.5)
hist_mids <- matrix(nrow = n_draw, ncol = max_breaks*1.5)
block_maxima <- matrix(nrow = n_draw, ncol = n_preds)
area_sums <- matrix(nrow = n_draw, ncol = n_preds)

n_events <- matrix(nrow = n_draw, ncol = n_preds)
for (i in 1:n_draw) {
  is_nb <- rbinom(n_preds, size = 1, prob = 1 - plogis(post$mu_full[i, 5, ]))
  n_events[i, ] <- is_nb * rnbinom(n_preds, size = exp(post$mu_full[i, 3, ]), mu = exp(post$mu_full[i, 4, ]))
}

assert_that(sum(is.na(n_events)) == 0)
gc()

min_fire_size <- rep(NA, n_draw)
burn_area_list <- list()

pb <- txtProgressBar(max = n_draw, style = 3)
for (i in 1:n_draw) {
  ba_vec <- list() # to hold all fire sizes for iteration i
  counter <- 1
  for (j in 1:n_preds) {
    if (n_events[i, j] > 0) {
      burn_areas <- 1e3 +
        weibull_scale_adj * rlomax(n = n_events[i, j],
                                   lambda = 1 / exp(post$mu_full[i, 1, j]),
                                   kappa = exp(post$mu_full[i, 2, j]))
      block_maxima[i, j] <- max(burn_areas)
      area_sums[i, j] <- sum(burn_areas)
      ba_vec[[counter]] <- burn_areas
      counter <- counter + 1
    }
  }
  ba_vec <- unlist(ba_vec)
  burn_area_list[[i]] <- ba_vec
  min_fire_size[i] <- min(ba_vec)
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
rm(ppred_df)
gc()

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

ba_df <-map(burn_area_list, as_tibble) %>%
  bind_rows(.id = 'iteration')

write_rds(ba_df, 'ppc_ba.rds')
