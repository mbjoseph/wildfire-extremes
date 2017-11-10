
# Create tables for LOO-IC for each fit -----------------------------------
library(loo)
library(tidyverse)
library(cowplot)

# possibly fetch from Amazon
system("aws s3 cp s3://earthlab-mjoseph .. --recursive --exclude '*' --include '*.rds'")

model_fits <- list.files(pattern = '*.fit.*\\.rds')

loo_df <- list()

for (i in seq_along(model_fits)) {
  fit_obj <- readRDS(model_fits[i])
  loglik_burn_area <- extract_log_lik(fit_obj, parameter_name = 'loglik_f')
  loglik_counts <- extract_log_lik(fit_obj, parameter_name = 'loglik_c')
  combined_loglik <- cbind(loglik_burn_area, loglik_counts)
  loo_total <- loo(combined_loglik, cores = 2)

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
    select(model, looic, se_looic)

  # clean up
  rm(fit_obj, loglik_burn_area, loo_total, loglik_counts)
  gc()
}

bind_rows(loo_df) %>%
  arrange(looic) %>%
  write_csv('loo-ic-table.csv')
