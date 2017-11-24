
# Create tables for LOO-IC for each fit -----------------------------------
library(loo)
library(tidyverse)
library(cowplot)
library(extraDistr)

source('R/02-explore.R')
source('R/make-stan-d.R')

# possibly fetch from Amazon
#system("aws s3 cp s3://earthlab-mjoseph .. --recursive --exclude '*' --include '*.rds'")

model_fits <- list.files(pattern = '*.fit.*\\.rds')

loo_df <- list()

holdout_counts <- count_df %>%
  mutate(er_ym = paste(NA_L3NAME, ym, sep = "_")) %>%
  filter(year >= cutoff_year)

holdout_indices <- match(holdout_counts$er_ym, st_covs$er_ym)

holdout_burns <-  mtbs %>%
  filter(FIRE_YEAR >= cutoff_year) %>%
  left_join(st_covs)

burn_indices <- match(holdout_burns$er_ym, holdout_counts$er_ym)

loglik_c <- list()
loglik_b <- list()

for (i in seq_along(model_fits)) {
  fit_obj <- readRDS(model_fits[i])
  loglik_burn_area <- extract_log_lik(fit_obj, parameter_name = 'loglik_f')
  loglik_counts <- extract_log_lik(fit_obj, parameter_name = 'loglik_c')
  combined_loglik <- cbind(loglik_burn_area, loglik_counts)
  loo_total <- loo(combined_loglik, cores = 1)

  loo_table_name <- gsub(pattern = 'fit', replacement = 'loo', x = model_fits[i])
  loo_table_name <- gsub(pattern = 'rds', replacement = 'csv', x = loo_table_name)
  loo_plot_name <- gsub(pattern = 'csv', replacement = 'png', x = loo_table_name)

  # save out a plot of each loo object
  model_type <- case_when(grepl('^gfit', model_fits[i]) ~ 'Gamma',
                          grepl('^wfit', model_fits[i]) ~ 'Weibull',
                          grepl('^lnfit', model_fits[i]) ~ 'Log-Normal',
                          grepl('^gpdfit', model_fits[i]) ~ 'Generalized Pareto',
                          grepl('^zinbgpdfit', model_fits[i]) ~ 'ZINB Generalized Pareto')
  plot_title <- paste('PSIS shape parameters:', model_type, 'model')
  png(filename = loo_plot_name, width = 12, height = 8, units = 'in', res = 300)
  plot(loo_total)
  title(plot_title)
  dev.off()
  # save a table of the loo values
  elements_to_save <- c('looic', 'se_looic')
  loo_df[[i]] <- loo_total[elements_to_save] %>%
    as_tibble %>%
    mutate(model = model_type) %>%
    dplyr::select(model, looic, se_looic) %>%
    mutate(i = i)

  # compute out of sample log likelihood -------------------------------------------------
  mu_full <- rstan::extract(fit_obj, pars = 'mu_full')[[1]]
  mu_holdout <- mu_full[, , holdout_indices]
  ll_c <- rep(NA, dim(mu_holdout)[1])
  ll_b <- rep(NA, dim(mu_holdout)[1])

  for (j in seq_along(ll_c)) {
    if (grepl('ZINB', model_type)) {
      # zinb count
      ll_c[j] <- ifelse(holdout_counts$n_fire == 0,
                        log(plogis(c(mu_holdout[j, 5, ])) +
                                  (1 - plogis(c(mu_holdout[j, 5, ]))) *
                          dnbinom(holdout_counts$n_fire,
                                  size = exp(c(mu_holdout[j, 3, ])),
                                  mu = exp(c(mu_holdout[j, 4, ])))),
                        # if not 0, then must be from ninom
                        log((1 - plogis(c(mu_holdout[j, 5, ])))) +
                          dnbinom(holdout_counts$n_fire,
                                size = exp(c(mu_holdout[j, 3, ])),
                                mu = exp(c(mu_holdout[j, 4, ])),
                                log = TRUE)) %>%
        mean()
    } else {
      # nb count
      ll_c[j] <- dnbinom(holdout_counts$n_fire,
                                     size = exp(c(mu_holdout[j, 3, ])),
                                     mu = exp(c(mu_holdout[j, 4, ])),
                         log = TRUE) %>%
        mean()
    }
    # compute log likelihood for heldout burn areas
    holdout_y <- (holdout_burns$R_ACRES - 1e3) / weibull_scale_adj
    if (grepl('Gamma', model_type)) {
      ll_b[j] <- dgamma(holdout_y,
                        shape = exp(c(mu_holdout[j, 1, burn_indices])),
                        rate = exp(c(mu_holdout[j, 2, burn_indices])),
                        log = TRUE) %>%
        mean()
    } else if (grepl('Weibull', model_type)) {
      ll_b[j] <- dweibull(holdout_y,
                        shape = exp(c(mu_holdout[j, 1, burn_indices])),
                        scale = exp(c(mu_holdout[j, 2, burn_indices])),
                        log = TRUE) %>%
        mean()
    } else if (grepl('Log-Normal', model_type)) {
      ll_b[j] <- dlnorm(holdout_y,
                        meanlog = exp(c(mu_holdout[j, 1, burn_indices])),
                        sdlog = exp(c(mu_holdout[j, 2, burn_indices])),
                        log = TRUE) %>%
        mean()
    } else if (grepl('Generalized Pareto', model_type)) {
      ll_b[j] <- dlomax(holdout_y,
                        lambda = 1 / exp(c(mu_holdout[j, 1, burn_indices])),
                        kappa = exp(c(mu_holdout[j, 2, burn_indices])),
                        log = TRUE) %>%
        mean()
    }
  }

  loglik_c[[i]] <- ll_c
  loglik_b[[i]] <- ll_b

  # clean up
  rm(fit_obj, loglik_burn_area, loo_total, loglik_counts, mu_full, mu_holdout)
  gc()
}

bind_rows(loo_df) %>%
  arrange(looic) %>%
  write_csv('loo-ic-table.csv')

model_names <- lapply(loo_df, function(x) paste(x$model, x$i)) %>%
  unlist()

names(loglik_c) <- model_names
plot(density(loglik_c[[1]]), xlim = c(-.56, -.4), ylim = c(0, 60))
for (i in 2:length(loglik_c)) {
  lines(density(loglik_c[[i]]), col = i)
}
legend('topleft', legend = model_names, lty = 1, col = seq_along(loglik_c), bty = 'n')

names(loglik_b) <- model_names
plot(density(log(-loglik_b[[1]])), xlim = c(-1, 0), ylim = c(0, 30))
for (i in 2:length(loglik_b)) {
  lines(density(log(-loglik_b[[i]]), na.rm = TRUE), col = i)
}
legend('topright', legend = model_names, lty = 1, col = seq_along(loglik_c), bty = 'n')
