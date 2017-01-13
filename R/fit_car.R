# CAR models --------------------------------------------------------------
library(rstan)
library(shinystan)
library(lme4)
source('R/explore.R')

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
         reg = as.numeric(factor(us_l3name)),
         year = fire_year + 1 - min(fire_year))

data_summary %>%
  ggplot(aes(x = year, y = n_fires / area, group = us_l3name)) +
  geom_line(alpha = .3) +
  scale_y_log10()

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
               year = data_summary$year)
m_init <- stan_model("stan/spatiotemporal-scale-free.stan", auto_write = TRUE)
m_fit <- sampling(m_init, data = stan_d, cores = 2, chains = 2,
                  iter = 300, init_r = .1)
traceplot(m_fit, pars = c("sigma_phi", "mu_sd", "sd_sd"))
traceplot(m_fit, pars = "gamma")
plot(m_fit, pars = "phi")
plot(m_fit, pars = "phiR")


# visualize spatiotemporal random effects
st <- m_fit %>%
  rstan::extract() %>%
  `[[`("phi") %>%
  reshape2::melt(varnames = c("iter", "reg", "year")) %>%
  group_by(year, reg) %>%
  summarize(mean = mean(value))

st %>%
  ggplot(aes(x = year, y = exp(mean), group = reg)) +
  geom_line() +
  scale_y_log10()

st %>%
  full_join(data_summary) %>%
  ungroup() %>%
  ggplot(aes(x = exp(mean), y = n_fires / area)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")

st %>%
  full_join(data_summary) %>%
  ungroup() %>%
  dplyr::select(us_l3name, fire_year, mean) %>%
  full_join(er_df) %>%
  dplyr::select(long, lat, group, mean, fire_year) %>%
  ggplot(aes(long, lat, group = group)) +
  geom_polygon(aes(fill = mean)) +
  facet_wrap(~ fire_year) +
  coord_equal() +
  scale_fill_viridis() +
  theme_map
