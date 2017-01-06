# CAR models --------------------------------------------------------------
library(rstan)
library(shinystan)
library(lme4)
source('R/explore.R')

# 1. Simple model for counts in each region -------------------------------
data_summary <- d %>%
  group_by(us_l3name) %>%
  summarize(n_fires = n()) %>%
  left_join(a_df) %>%
  full_join(all_ers) %>%
  group_by(us_l3name) %>%
  summarize(n_fires = unique(n_fires),
            area = unique(area)) %>%
  mutate(n_fires = ifelse(is.na(n_fires), 0, n_fires))

stopifnot(all(rownames(W) == data_summary$us_l3name))

stan_d <- list(n = nrow(data_summary),
               p = 1,
               X = matrix(1, nrow(data_summary), ncol = 1),
               y = data_summary$n_fires,
               log_offset = log(data_summary$area),
               W_n = sum(W) / 2,
               W = W)

m_init <- stan_model("stan/sparse_car.stan", auto_write = TRUE)
m_fit <- sampling(m_init, data = stan_d, cores = 2, chains = 2, iter = 2000)
m_fit
traceplot(m_fit)

data_summary$ranef <- m_fit %>%
  rstan::extract() %>%
  `[[`('phi') %>%
  apply(2, median)

data_summary %>%
  dplyr::select(us_l3name, ranef) %>%
  full_join(er_df) %>%
  ggplot(aes(long, lat, group = group)) +
  geom_polygon(aes(fill = ranef)) +
  coord_equal() +
  scale_fill_viridis() +
  theme_map


# 2. Regionally varying slopes and intercepts -----------------------------
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

stan_d <- list(n = nrow(data_summary),
               p = 1,
               X = matrix(1, nrow(data_summary)),
               y = data_summary$n_fires,
               J = max(data_summary$reg),
               reg = data_summary$reg,
               t = data_summary$cyear,
               log_offset = log(data_summary$area),
               W_n = sum(W) / 2,
               W = W,
               n_year = length(unique(data_summary$fire_year)),
               year = data_summary$year)

m_init <- stan_model("stan/2-varying-slope.stan", auto_write = TRUE)
m_fit <- sampling(m_init, data = stan_d, cores = 2, chains = 2, iter = 3000)
m_fit
traceplot(m_fit, pars = c("beta", "tau_phi", 'tau_phi1', "gamma"))
plot(m_fit, pars = "phi")




st <- m_fit %>%
  rstan::extract() %>%
  `[[`("phi") %>%
  reshape2::melt(varnames = c("iter", "year", "reg")) %>%
  group_by(year, reg) %>%
  summarize(mean = mean(value))

st %>%
  ggplot(aes(x = year, y = mean, group = reg)) +
  geom_line()



st %>%
  full_join(data_summary) %>%
  ungroup() %>%
  dplyr::select(us_l3name, fire_year, mean) %>%
  full_join(er_df) %>%
  ggplot(aes(long, lat, group = group)) +
  geom_polygon(aes(fill = mean)) +
  facet_wrap(~ fire_year) +
  coord_equal() +
  scale_fill_viridis() +
  theme_map


