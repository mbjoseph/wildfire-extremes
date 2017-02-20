# CAR models --------------------------------------------------------------
library(rstan)
library(plotly)
library(scatterD3)

source('R/explore.R')

# get count data
data_summary <- d %>%
  group_by(us_l3name, fire_year) %>%
  summarize(n_fires = n()) %>%
  left_join(a_df) %>%
  full_join(all_ers) %>%
  ungroup() %>%
  dplyr::select(-area, -n_neighbors) %>%
  complete(us_l3name, fire_year,
           fill = list(n_fires = 0)) %>%
  full_join(all_ers) %>%
  mutate(cyear = c(scale(fire_year)),
         year = fire_year + 1 - min(fire_year),
         freg = factor(us_l3name,
                       levels = levels(factor(all_ers$us_l3name))),
         reg = as.numeric(freg))

# get fire size data
fire_sizes <- d %>%
  dplyr::select(us_l3name, fire_year, fire_size) %>%
  mutate(cyear = c(scale(fire_year)),
         freg = factor(us_l3name,
                                 levels = levels(data_summary$freg)),
         reg = as.numeric(freg),
         year = fire_year + 1 - min(fire_year)) %>%
  arrange(us_l3name, fire_year)



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

if (fit_car) {
  m_init <- stan_model("stan/spatiotemporal-scale-free.stan", auto_write = TRUE)

  m_fit <- sampling(m_init, data = stan_d, cores = 4, chains = 4,
                    iter = 300, init_r = .0001,
                    pars = c("sigma_phi_z", "sigma_phi_y", "sd_sd",
                             "gamma",
                             "phi_z", "phi_y", "phi_zR", "phi_yR", "sigma_z",
                             "z_j", "z0", "sigma_z_j",
                             "alpha",
                             "beta_phi_y", "mu_beta_phi", "sigma_beta_phi",
                             "zmax_new", "y_new", "log_lambda_new", "mu_z_new",
                             "area_burned"))
  saveRDS(m_fit, "m_fit.rds")
}
