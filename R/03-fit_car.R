# CAR models --------------------------------------------------------------
library(rstan)
library(plotly)
library(scatterD3)

source('R/01-clean_data.R')
source('R/02-explore.R')


# bundle data for stan
stan_d <- list(n = nrow(data_summary),
               x = data_summary$cyear,
               y = data_summary$n_fires,
               J = max(data_summary$reg),
               reg = data_summary$reg,
               t = data_summary$cyear,
               log_offset = log(all_ers$area),
               W_n = sum(W) / 2,
               W = W,
               n_year = length(unique(data_summary$year)) + 3,
               year = data_summary$year,
               nz = nrow(fire_sizes),
               z = log(fire_sizes$fire_size),
               reg_z = fire_sizes$reg,
               year_z = fire_sizes$year)

# fit model
m_init <- stan_model("stan/spatiotemporal-scale-free.stan", auto_write = TRUE)

m_fit <- sampling(m_init, data = stan_d, cores = 3, chains = 3,
                  iter = 400, init_r = .0001,
                  pars = c("sigma_phi_z", "sigma_phi_y", "sd_sd",
                           "gamma",
                           "phi_z", "phi_y", "phi_zR", "phi_yR", "sigma_z",
                           "z0",
                           "alpha",
                           "beta_phi_y", "mu_beta_phi", "sigma_beta_phi",
                           "zmax_new", "y_new", "log_lambda_new", "mu_z_new",
                           "area_burned"))
saveRDS(m_fit, "m_fit.rds")
