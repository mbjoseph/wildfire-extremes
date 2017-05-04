# CAR models --------------------------------------------------------------
source('R/02-explore.R')

library(rstan)
library(plotly)
library(scatterD3)

# bundle data for stan
stan_d <- list(n = nrow(data_summary),
               y = data_summary$n_fires,
               J = max(data_summary$reg),
               reg = data_summary$reg,
               t = data_summary$cyear,
               log_offset = log(all_ers$area / 1E8),
               W_n = sum(W) / 2,
               W = W,
               n_year = length(unique(data_summary$num_ym)),
               year = data_summary$num_ym,
               nz = nrow(fire_sizes),
               z = log(fire_sizes$fire_size),
               reg_z = fire_sizes$reg,
               year_z = fire_sizes$num_ym)

# fit model
m_init <- stan_model("stan/spatiotemporal-scale-free.stan", auto_write = TRUE)

m_fit <- sampling(m_init, data = stan_d, cores = 2, chains = 2,
                  iter = 10, init_r = .1,
                  pars = c("sigma_phi_z", "sigma_phi_y", "sd_sd",
                           "gamma",
                           "phi_z", "phi_y", "phi_zR", "phi_yR", "sigma_z",
                           "z0",
                           "alpha",
                           "beta_phi_y",
                           "zmax_new", "y_new", "log_lambda_new", "mu_z_new",
                           "area_burned"
                           ))
saveRDS(m_fit, "m_fit.rds")


traceplot(m_fit, pars = c("sigma_z", 'alpha'))
