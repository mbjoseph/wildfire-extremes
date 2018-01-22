
# Count plots -------------------------------------------------------------

source('R/02-explore.R')
source('R/make-stan-d.R')
library(ggridges)
library(viridis)
library(ggthemes)
library(plotly)

fit <- read_rds(path = list.files(pattern = 'zinb_fit.*'))

# Evaluate convergence ----------------------------------------------------
traceplot(fit, inc_warmup = TRUE)

traceplot(fit, pars = c('tau', 'c', 'alpha', 'Rho_beta'))

# Extract posterior draws and visualize results ---------------------------
post <- rstan::extract(fit)
str(post)
rm(fit)
gc()



# Coefficients that seem to be far from zero ------------------------------
beta_df <- post$beta %>%
  reshape2::melt(varnames = c('iter', 'dim', 'col')) %>%
  tbl_df

beta_summary <- beta_df %>%
  group_by(dim, col) %>%
  summarize(median = median(value),
            lo = quantile(value, 0.025),
            hi = quantile(value, 0.975),
            p_neg = mean(value < 0),
            p_pos = mean(value > 0)) %>%
  ungroup %>%
  mutate(variable = colnamesX[col],
         nonzero = p_neg > .9 | p_pos > .9)

beta_summary %>%
  filter(nonzero) %>%
  select(dim, col, p_neg, p_pos, variable) %>%
  left_join(beta_df) %>%
  ggplot(aes(value, variable)) +
  geom_density_ridges() +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  facet_wrap(~ dim, scales = 'free') +
  theme(axis.text.y = element_text(size = 7)) +
  ggtitle('Estimated effects on # of fires')
ggsave('fig/fire-effs.png', width = 8, height = 8)

beta_summary %>%
  filter(nonzero) %>%
  ggplot(aes(x = median, y = reorder(variable, median))) +
  geom_point() +
  facet_wrap(~dim, scales = 'free') +
  geom_errorbarh(aes(xmin = lo, xmax = hi)) +
  theme(axis.text.y = element_text(size = 8))
ggsave('fig/count-effs.png', width = 6, height = 12)

rm(beta_df)
gc()

mu_df <- post$mu_full %>%
  reshape2::melt(varnames = c('iter', 'dim', 'id')) %>%
  as_tibble() %>%
  group_by(id, dim) %>%
  summarize(med = median(value)) %>%
  spread(dim, med) %>%
  mutate(expected_val = plogis(`2`) * exp(`1`))

st_covs <- st_covs %>%
  group_by(NA_L3NAME) %>%
  mutate(mean_hd = mean(housing_density))

st_covs %>%
  full_join(mu_df) %>%
  ggplot(aes(x = rmin, y = expected_val / (1000 * area), color = housing_density)) +
  geom_point(size = .5) +
  theme_minimal() +
  facet_wrap(~ NA_L2NAME,
             labeller = labeller(NA_L2NAME = label_wrap_gen(16))) +
  theme(strip.text = element_text(size = 8)) +
  scale_color_viridis(trans = 'log', 'Housing density') +
  scale_y_log10() +
  xlab('Mean daily minimum humidity') +
  ylab('Expected fire density (# per sq. km)')

st_covs %>%
  full_join(mu_df) %>%
  ggplot(aes(x = vs, y = expected_val / (1000 * area), color = housing_density)) +
  geom_point(size = .5, alpha = .5) +
  theme_minimal() +
  facet_wrap(~ NA_L2NAME,
             labeller = labeller(NA_L2NAME = label_wrap_gen(16))) +
  theme(strip.text = element_text(size = 8)) +
  scale_color_viridis(trans = 'log', 'Housing density') +
  xlab('Mean daily wind speed') +
  ylab('Expected fire density (# per sq. km)') +
  scale_y_log10()

st_covs %>%
  full_join(mu_df) %>%
  ggplot(aes(x = prev_12mo_precip, y = expected_val / (1000 * area), color = housing_density)) +
  geom_point(size = .5, alpha = .5) +
  theme_minimal() +
  facet_wrap(~ NA_L2NAME,
             labeller = labeller(NA_L2NAME = label_wrap_gen(16))) +
  theme(strip.text = element_text(size = 8)) +
  scale_color_viridis(trans = 'log', 'Housing density') +
  xlab('Mean previous 12 month precipitation') +
  ylab('Expected fire density (# per sq. km)') +
  scale_x_log10() +
  scale_y_log10()

st_covs %>%
  full_join(mu_df) %>%
  ggplot(aes(x = pr, y = expected_val / (1000 * area), color = housing_density)) +
  geom_point(size = .5, alpha = .5) +
  theme_minimal() +
  facet_wrap(~ NA_L2NAME,
             labeller = labeller(NA_L2NAME = label_wrap_gen(16))) +
  theme(strip.text = element_text(size = 8)) +
  scale_color_viridis(trans = 'log', 'Housing density') +
  xlab('Same month precipitation') +
  ylab('Expected fire density (# per sq. km)') +
  scale_y_log10()

st_covs %>%
  full_join(mu_df) %>%
  ggplot(aes(x = tmmx, y = expected_val / (1000 * area), color = housing_density)) +
  geom_point(size = .5, alpha = .5) +
  theme_minimal() +
  facet_wrap(~ NA_L2NAME,
             labeller = labeller(NA_L2NAME = label_wrap_gen(16))) +
  theme(strip.text = element_text(size = 8)) +
  scale_color_viridis(trans = 'log', 'Housing density') +
  xlab('Air temperature') +
  ylab('Expected fire density (# per sq. km)') +
  scale_y_log10()


# Predicting the expected value across a range of vals --------------------
ng <- 100
st_pred <- expand.grid(ctri = 0,
                       cpr = 0,
                       ctmx = 0,
                       cpr12 = 0,
                       cvs = 0,
                       chd = seq(min(st_covs$chd), max(st_covs$chd), length.out = ng),
                       crmin = seq(min(st_covs$crmin), max(st_covs$crmin), length.out = ng),
                       NA_L3NAME = c('Acadian Plains and Hills', 'dummy'),
                       NA_L2NAME = c('MIXED WOOD PLAINS', 'dummy'),
                       NA_L1NAME = c('EASTERN TEMPERATE FORESTS', 'dummy')) %>%
  as_tibble

range(st_covs$chd)
range(st_pred$chd)

X_bs_pred <- list()
X_bs_pred_df <- list()
for (i in seq_along(vars)) {
  X_bs_pred[[i]] <- predict(X_bs[[i]], newx = st_pred[[vars[i]]])
  X_bs_pred_df[[i]] <- X_bs_pred[[i]] %>%
    as_tibble
  names(X_bs_pred_df[[i]]) <- paste('bs', vars[[i]], 1:df_each, sep = '_')
}
X_bs_pred_df <- bind_cols(X_bs_pred_df)
assert_that(!any(is.na(X_bs_pred_df)))

# merge back into prediction df
st_pred <- st_pred %>%
  bind_cols(lapply(X_bs_pred_df, c)) %>%
  as_tibble

plot(st_pred$crmin, st_pred$bs_crmin_1)
points(st_covs$crmin, st_covs$bs_crmin_1,
       col = scales::alpha(2, .1), cex = .5)

plot(st_pred$crmin, st_pred$bs_crmin_4)
points(st_covs$crmin, st_covs$bs_crmin_4,
       col = scales::alpha(2, .1), cex = .5)

plot(st_pred$cpr, st_pred$bs_cpr_4, ylim = c(0, 1), xlim = c(-4, 4))
points(st_covs$cpr, st_covs$bs_cpr_4,
       col = scales::alpha(2, .1), cex = .5)


# construct the predictive design matrix
X_pred <- model.matrix(as.formula(paste('~ 0 + ctri + ', l1_terms, l2_terms, l3_terms, sep = ' + ')),
                    data = st_pred)

# prune dummy levels that were necessary to get contrasts
X_pred <- X_pred[!(st_pred$NA_L1NAME == 'dummy' | st_pred$NA_L2NAME == 'dummy' | st_pred$NA_L3NAME == 'dummy'), ]
X_pred <- X_pred[, !grepl('dummy',colnames(X_pred))]

# which columns to pull from?
x_idx <- match(colnames(X_pred), colnames(X)) %>%
  na.omit %>%
  c
# add predicted log_lambda as new column (use average area)
pred_mat <- array(dim = c(nrow(st_pred), 2, length(post$lp__)))
for (i in 1:2) {
  for (j in 1:length(post$lp__)) {
    pred_mat[, i, j] <- X_pred %*% post$beta[j, i, x_idx] + mean(stan_d$log_area)
  }
  if (i == 1) {
    pred_mat[, i, ] <- pred_mat[, i, ] + mean(stan_d$log_area)
  }
}

pred_summ <- array(dim = c(nrow(st_pred), 2, 3))
quantiles <- c(.05, .5, .95)
qnames <- c('lo', 'med', 'hi')
for (i in 1:nrow(st_pred)) {
  for (j in 1:2) {
    pred_summ[i, j, ] <- quantile(pred_mat[i, j, ], quantiles)
  }
}

pred_df <- reshape2::melt(pred_summ, varnames = c('pred_idx', 'dim', 'q')) %>%
  as_tibble %>%
  mutate(quantile = qnames[q]) %>%
  dplyr::select(-q)

st_pred$pred_idx <- 1:nrow(st_pred)

all_preds <- st_pred %>%
  dplyr::select(pred_idx, chd, crmin) %>%
  full_join(pred_df)

# expected value from neg binom
all_preds %>%
  filter(dim == 1) %>%
  spread(quantile, value) %>%
  ggplot(aes(crmin, exp(med), color = chd, group = chd, fill = chd)) +
  theme_minimal() +
  scale_color_viridis_c() +
  geom_line() +
  scale_fill_viridis_c()


# expected value from neg binom
all_preds %>%
  filter(dim == 2) %>%
  spread(quantile, value) %>%
  ggplot(aes(crmin, plogis(med), color = chd, group = chd, fill = chd)) +
  theme_minimal() +
  scale_color_viridis_c() +
  geom_line() +
  scale_fill_viridis_c()




plot_ly(filter(all_preds, dim == 1, quantile == 'med'),
        x = ~crmin, y = ~chd,
        z = ~exp(value),
        color = ~exp(value),
        marker = list(size = 2)) %>%
  add_markers()


plot_ly(filter(all_preds, dim == 2, quantile == 'med'),
        x = ~crmin, y = ~chd,
        z = ~plogis(value),
        color = ~plogis(value),
        marker = list(size = 2)) %>%
  add_markers()

plot_ly(filter(all_preds, dim == 2, quantile == 'med'),
        x = ~crmin, y = ~exp(chd * sd(log(st_covs$housing_density)) +
                               mean(log(st_covs$housing_density))),
        z = ~plogis(value)) %>%
  add_trace(type="mesh3d") %>%
  add_trace(type = 'mesh3d', data = filter(all_preds, dim == 2, quantile == 'lo')) %>%
  add_trace(type = 'mesh3d', data = filter(all_preds, dim == 2, quantile == 'hi'))


