# CAR models --------------------------------------------------------------
library(rstan)
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
               n_year = length(unique(data_summary$year)) + 5,
               year = data_summary$year,
               nz = nrow(fire_sizes),
               z = log(fire_sizes$fire_size),
               reg_z = fire_sizes$reg,
               year_z = fire_sizes$year)

m_init <- stan_model("stan/spatiotemporal-scale-free.stan", auto_write = TRUE)

m_fit <- sampling(m_init, data = stan_d, cores = 2, chains = 2,
                  iter = 400, init_r = .001,
                  pars = c("sigma_phi_z", "sigma_phi_y", "mu_sd", "sd_sd",
                           "gamma",
                           "phi_z", "phi_y", "phi_zR", "phi_yR", "sigma_z",
                           "z_j", "z0", "sigma_z_j",
                           "alpha",
                           "beta_phi_y", "mu_beta_phi", "sigma_beta_phi",
                           "zmax_new", "y_new", "log_lambda_new", "mu_z_new",
                           "area_burned"))


# Evaluate convergence ----------------------------------------------------

traceplot(m_fit, pars = "sigma_phi_z")
traceplot(m_fit, pars = "sigma_phi_y")
traceplot(m_fit, pars = c("mu_sd", "sd_sd"))
traceplot(m_fit, pars = "gamma")
traceplot(m_fit, pars = c("z0", "sigma_z_j"))
traceplot(m_fit, pars = "alpha")
traceplot(m_fit, pars = "sigma_z")
traceplot(m_fit, pars = c("mu_beta_phi", "sigma_beta_phi"))



plot(m_fit, pars = "phi_zR") +
  coord_flip()


plot(m_fit, pars = "phi_yR") +
  coord_flip()

plot(m_fit, pars = "phi_y") +
  coord_flip()

plot(m_fit, pars = "beta_phi_y")


# compare observed numbers to predicted
m_fit %>%
  rstan::extract() %>%
  `[[`("y_new") %>%
  reshape2::melt(varnames = c("iter", "reg", "year")) %>%
  group_by(year, reg) %>%
  summarize(median = median(value),
            lo = quantile(value, 0.025),
            hi = quantile(value, 0.975)) %>%
  mutate(par = "y_new") %>%
  full_join(data_summary) %>%
  ggplot(aes(x = fire_year)) +
  geom_ribbon(aes(ymin = lo, ymax = hi),
              alpha = .5, color = NA, fill = "dodgerblue") +
  geom_line(aes(y = median), color = "dodgerblue") +
  geom_point(aes(y = n_fires)) +
  facet_wrap(~ us_l3name) +
  scale_y_log10() +
  theme_minimal()



# Posterior predictive check for maxima -----------------------------------
maxima <- fire_sizes %>%
  group_by(us_l3name, fire_year) %>%
  summarize(zmax = max(log(fire_size))) %>%
  full_join(data_summary)

ppc_max <- m_fit %>%
  rstan::extract() %>%
  `[[`("zmax_new") %>%
  reshape2::melt(varnames = c("iter", "reg", "year")) %>%
  filter(value > 0) %>%
  group_by(year, reg) %>%
  summarize(median = median(value),
            lo = quantile(value, 0.025),
            hi = quantile(value, 0.975)) %>%
  full_join(maxima) %>%
  mutate(Projected = year > length(unique(data_summary$year)))

ppc_max %>%
  ggplot(aes(x = year, group = reg)) +
  geom_ribbon(aes(ymin = exp(lo), ymax = exp(hi)),
              alpha = .05, fill = "dodgerblue") +
  geom_line(aes(y = exp(median)), alpha = .4) +
  theme_minimal()

ppc_max %>%
  ggplot(aes(x = year + min(data_summary$fire_year) - 1)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = Projected),
              alpha = .5) +
  geom_line(aes(y = median, color = Projected)) +
  geom_point(aes(y = zmax), size = .5) +
  theme_minimal() +
  facet_wrap(~ reg, scales = "free")


# Posterior predictive check for burn area --------------------------------
burn_areas <- fire_sizes %>%
  group_by(us_l3name, fire_year) %>%
  summarize(zsum = sum(log(fire_size))) %>%
  full_join(data_summary)

ppc_sum <- m_fit %>%
  rstan::extract() %>%
  `[[`("area_burned") %>%
  reshape2::melt(varnames = c("iter", "reg", "year")) %>%
  filter(value > 0) %>%
  group_by(year, reg) %>%
  summarize(median = median(value),
            lo = quantile(value, 0.025),
            hi = quantile(value, 0.975)) %>%
  full_join(burn_areas) %>%
  mutate(Projected = year > length(unique(data_summary$year)))

ppc_sum %>%
  ggplot(aes(x = year, group = reg)) +
  geom_ribbon(aes(ymin = lo, ymax = hi),
              alpha = .05, fill = "dodgerblue") +
  geom_line(aes(y = median), alpha = .4) +
  theme_minimal() +
  scale_y_log10()


ppc_sum %>%
  filter(!Projected) %>%
  ggplot(aes(x = year + min(data_summary$fire_year) - 1)) +
  geom_ribbon(aes(ymin = lo, ymax = hi),
              alpha = .5, fill = "dodgerblue") +
  geom_line(aes(y = median)) +
  geom_point(aes(y = zsum), size = .5) +
  theme_minimal() +
  facet_wrap(~ reg, scales = "free")



st_y <- m_fit %>%
  rstan::extract() %>%
  `[[`("log_lambda_new") %>%
  reshape2::melt(varnames = c("iter", "reg", "year")) %>%
  group_by(year, reg) %>%
  summarize(median = median(value),
            lo = quantile(value, 0.025),
            hi = quantile(value, 0.975)) %>%
  mutate(par = "log_lambda")


st_y %>%
  ggplot(aes(x = year, y = median, group = reg)) +
  geom_line()


st_z <- m_fit %>%
  rstan::extract() %>%
  `[[`("mu_z_new") %>%
  reshape2::melt(varnames = c("iter", "reg", "year")) %>%
  group_by(year, reg) %>%
  summarize(median = median(value),
            lo = quantile(value, 0.025),
            hi = quantile(value, 0.975)) %>%
  mutate(par = "mu_z")

st_z %>%
  ggplot(aes(x = year, y = median, group = reg)) +
  geom_line()


full_join(st_y, st_z) %>%
  dplyr::select(-lo, -hi) %>%
  spread(par, median) %>%
  ggplot(aes(x = log_lambda, y = mu_z, group = reg, color = year)) +
  geom_path() +
  scale_color_viridis() +
  theme_minimal() +
  facet_wrap(~ reg)




# visualize spatiotemporal random effects
st_y <- m_fit %>%
  rstan::extract() %>%
  `[[`("phi_y") %>%
  reshape2::melt(varnames = c("iter", "reg", "year")) %>%
  group_by(year, reg) %>%
  summarize(median = median(value),
            lo = quantile(value, 0.025),
            hi = quantile(value, 0.975)) %>%
  mutate(par = "y")


st_y %>%
  ggplot(aes(x = year, y = exp(median), group = reg)) +
  geom_line() +
  scale_y_log10()

st_y %>%
  full_join(data_summary) %>%
  ungroup() %>%
  ggplot(aes(x = exp(median), y = n_fires / (area))) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")


st_z <- m_fit %>%
  rstan::extract() %>%
  `[[`("phi_z") %>%
  reshape2::melt(varnames = c("iter", "reg", "year")) %>%
  group_by(year, reg) %>%
  summarize(median = median(value),
            lo = quantile(value, 0.025),
            hi = quantile(value, 0.975)) %>%
  mutate(par = "z")


st_z %>%
  ggplot(aes(x = year, y = exp(median), group = reg)) +
  geom_line() +
  scale_y_log10() +
#  geom_ribbon(aes(ymin = exp(lo), ymax = exp(hi)), alpha = .03) +
  theme_minimal()



m_fit %>%
  rstan::extract() %>%
  `[[`("z_t") %>%
  reshape2::melt(varnames = c("iter", "year")) %>%
  group_by(year) %>%
  summarize(median = median(value),
            lo = quantile(value, 0.025),
            hi = quantile(value, 0.975)) %>%
  ggplot(aes(x = year + 1991, y = median)) +
  geom_line() +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .5) +
  theme_minimal() +
  xlab("Year") +
  ylab("Year-specific adjustment")


# Visualize link between fire occurrence and burn area
beta_phi_y <- m_fit %>%
  rstan::extract() %>%
  `[[`("beta_phi_y") %>%
  reshape2::melt(varnames = c("iter", "reg")) %>%
  group_by(reg) %>%
  summarize(median = median(value),
            lo = quantile(value, 0.025),
            hi = quantile(value, 0.975)) %>%
  mutate(over_zero = lo < 0 & hi > 0)

beta_phi_y %>%
  left_join(all_ers) %>%
  ggplot(aes(y = reorder(us_l3name, median), x = median, color = over_zero)) +
  geom_point() +
  geom_segment(aes(x = lo, xend = hi, yend = reorder(us_l3name, median))) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  ylab("Ecoregion") +
  xlab("Effect of fire occurence on burn area")

full_join(st_y, st_z) %>%
  ungroup() %>%
  dplyr::select(-lo, -hi) %>%
  spread(key = par, value = median) %>%
  left_join(all_ers) %>%
  ggplot(aes(x = y, y = z)) +
  xlab("Fire occurrence effects") +
  ylab("Fire burn area effects") +
  geom_point(alpha = .6) +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~us_l3name)


beta_phi_y %>%
  left_join(all_ers) %>%
  full_join(er_df) %>%
  ggplot(aes(x = long, y = lat, group = group)) +
  geom_polygon(aes(fill = median), color = NA) +
  coord_equal() +
  labs(x = "Longitude", y = "Latitude") +
  scale_fill_viridis("Effect of fire density on burn area") +
  theme_map

beta_phi_y %>%
  left_join(all_ers) %>%
  full_join(er_df) %>%
  ggplot(aes(x = long, y = lat, group = group)) +
  geom_polygon(aes(fill = hi - lo), color = NA) +
  coord_equal() +
  labs(x = "Longitude", y = "Latitude") +
  scale_fill_viridis("Uncertainty in effect of\nfire density on burn area") +
  theme_map
