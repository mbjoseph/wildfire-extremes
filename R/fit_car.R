# CAR models --------------------------------------------------------------
library(rstan)
library(shinystan)
library(lme4)
library(gganimate)
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
         freg = factor(us_l3name),
         reg = as.numeric(freg),
         year = fire_year + 1 - min(fire_year))

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
               log_offset = log(data_summary$area),
               W_n = sum(W) / 2,
               W = W,
               n_year = length(unique(data_summary$year)),
               year = data_summary$year,
               nz = nrow(fire_sizes),
               z = log(fire_sizes$fire_size),
               reg_z = fire_sizes$reg,
               year_z = fire_sizes$year)



m_init <- stan_model("stan/spatiotemporal-scale-free.stan", auto_write = TRUE)
m_fit <- sampling(m_init, data = stan_d, cores = 4, chains = 4,
                  iter = 2000, init_r = .01,
                  pars = c("sigma_phi_z", "sigma_phi_y", "mu_sd", "sd_sd", "gamma",
                           "phi_z", "phi_y", "phi_zR", "phi_yR", "sigma_z",
                           "z_t", "mu_z_t", "sigma_z_t", "alpha"))
traceplot(m_fit, pars = "sigma_phi_z")
traceplot(m_fit, pars = "sigma_phi_y")
traceplot(m_fit, pars = c("mu_sd", "sd_sd"))
traceplot(m_fit, pars = "gamma")
traceplot(m_fit, pars = c("mu_z_t", "sigma_z", "sigma_z_t"))
traceplot(m_fit, pars = "z_t")
traceplot(m_fit, pars = "alpha")

plot(m_fit, pars = "phi_zR") +
  coord_flip()


plot(m_fit, pars = "phi_yR") +
  coord_flip()



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
  scale_y_log10()


full_join(st_y, st_z) %>%
  ungroup() %>%
  dplyr::select(-lo, -hi) %>%
  spread(key = par, value = median) %>%
  ggplot(aes(x = y, y = z, color = factor(reg))) +
  xlab("Fire occurrence effects") +
  ylab("Fire burn area effects") +
  stat_smooth(method = "lm", alpha = .1) +
  theme_minimal() +
  theme(legend.position = "none")

# animated version
p <- full_join(st_y, st_z) %>%
  ungroup() %>%
  dplyr::select(-lo, -hi) %>%
  spread(key = par, value = median) %>%
  ggplot(aes(x = y, y = z, color = factor(reg), frame = year)) +
  xlab("Fire occurrence effects") +
  ylab("Fire burn area effects") +
  geom_point() +
  theme_minimal() +
  theme(legend.position = "none")
gganimate(p, filename = "out.gif")
