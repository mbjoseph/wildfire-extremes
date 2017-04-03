# CAR models --------------------------------------------------------------
library(rstan)
library(plotly)
library(scatterD3)

source('R/02-explore.R')

# bundle data for stan
stan_d <- list(n = nrow(data_summary),
               x = data_summary$cyear,
               y = data_summary$n_fires,
               J = max(data_summary$reg),
               reg = data_summary$reg,
               t = data_summary$cyear,
               log_offset = log(all_ers$area / 1E8),
               W_n = sum(W) / 2,
               W = W,
               n_year = length(unique(data_summary$year)),
               year = data_summary$year,
               nz = nrow(fire_sizes),
               z = log(fire_sizes$fire_size),
               prcp_z = log(fire_sizes$prcp) / 1.5,
               reg_z = fire_sizes$reg,
               year_z = fire_sizes$year,
               prcp_y = log(data_summary$prcp) / 1.5)

# fit model
m_init <- stan_model("stan/spatiotemporal-scale-free.stan", auto_write = TRUE)

m_fit <- sampling(m_init, data = stan_d, cores = 3, chains = 3,
                  iter = 800, init_r = .1,
                  pars = c("sigma_phi_z", "sigma_phi_y", "sd_sd",
                           "gamma", "sd_gamma", "mu_gamma",
                           "phi_z", "phi_y", "phi_zR", "phi_yR", "sigma_z",
                           "z0",
                           "alpha",
                           "beta_phi_y", "mu_beta_phi", "sigma_beta_phi",
                           "beta_prcp", "sigma_prcp", "mu_prcp"
                           # "zmax_new", "y_new", "log_lambda_new", "mu_z_new",
                           # "area_burned"
                           ))
saveRDS(m_fit, "m_fit.rds")


traceplot(m_fit, pars = c("mu_prcp", "sigma_prcp", 'alpha'))

prcp_eff <- rstan::extract(m_fit) %>%
  `[[`("beta_prcp") %>%
  reshape2::melt(varnames = c("iter", "par", "ecoregion")) %>%
  tbl_df %>%
  mutate(us_l3name = unique(data_summary$us_l3name)[ecoregion],
         effect = ifelse(par == 1, "Effect on fire density",
                         "Effect on mean burn area")) %>%
  group_by(us_l3name, effect) %>%
  summarize(med = median(value),
            lo = quantile(value, .1),
            hi = quantile(value, .9))

prcp_eff %>%
  mutate(nonzero = lo > 0 | hi < 0) %>%
  ggplot(aes(x = med, y = us_l3name, color = nonzero)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point() +
  facet_wrap(~effect) +
  geom_segment(aes(x = lo, xend = hi, yend = us_l3name)) +
  xlab("Estimate") +
  ylab("") +
  ggtitle("Estimated effects of precipitation on fire density and burn area") +
  theme_minimal()
