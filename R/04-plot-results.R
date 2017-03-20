library(rstan)
library(scatterD3)
source("R/01-clean_data.R")
source("R/02-explore.R")

m_fit <- readRDS("m_fit.rds")


min_year <- min(d$fire_year)
regions <- levels(data_summary$freg)

# Evaluate convergence ----------------------------------------------------

traceplot(m_fit, pars = "sigma_phi_z")
traceplot(m_fit, pars = "sigma_phi_y")
traceplot(m_fit, pars = c("sd_sd"))
traceplot(m_fit, pars = "gamma")
traceplot(m_fit, pars = c("z0"))
traceplot(m_fit, pars = "alpha")
traceplot(m_fit, pars = "sigma_z")
traceplot(m_fit, pars = c("mu_beta_phi", "sigma_beta_phi"))

plot(m_fit, pars = "phi_zR") +
  coord_flip()

plot(m_fit, pars = "phi_yR") +
  coord_flip()

plot(m_fit, pars = "phi_y") +
  coord_flip()

traceplot(m_fit, pars = "beta_phi_y")


# compare observed numbers to predicted
m_fit %>%
  rstan::extract() %>%
  `[[`("y_new") %>%
  reshape2::melt(varnames = c("iter", "reg", "year")) %>%
  group_by(year, reg) %>%
  summarize(median = median(value),
            lo = quantile(value, 0.025),
            hi = quantile(value, 0.975)) %>%
  mutate(par = "y_new",
         us_l3name = regions[reg],
         fire_year = year + min_year - 1) %>%
  mutate(Projected = year > length(unique(data_summary$year))) %>%
  ggplot(aes(x = fire_year, color = Projected, fill = Projected)) +
  geom_ribbon(aes(ymin = lo, ymax = hi),
              alpha = .5, color = NA) +
  geom_line(aes(y = median)) +
  geom_point(aes(x = fire_year, y = n_fires), data = data_summary,
             inherit.aes = FALSE) +
  facet_wrap(~ us_l3name) +
  scale_y_log10() +
  theme_minimal() +
  xlab("Year") +
  ylab("Number of fires")
ggsave("fig/ppc-n.pdf", width = 20, height = 10)



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
  ungroup() %>%
  mutate(year = year + min_year - 1,
         region = regions[reg]) %>%
  ggplot(aes(x = year)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = Projected),
              alpha = .5) +
  geom_line(aes(y = median, color = Projected)) +
  geom_point(aes(y = zmax), size = .5) +
  theme_minimal() +
  facet_wrap(~ region) +
  xlab("Year") +
  ylab("Maximum burn area")
ggsave("fig/ppc-max-burn-area.pdf", width = 20, height = 10)


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
  ungroup() %>%
  mutate(year = year + min_year - 1,
         region = regions[reg]) %>%
  ggplot(aes(x = year, fill = Projected)) +
  geom_ribbon(aes(ymin = lo, ymax = hi),
              alpha = .5) +
  geom_line(aes(y = median)) +
  geom_point(aes(y = zsum), size = .5) +
  theme_minimal() +
  facet_wrap(~ region) +
  scale_y_log10() +
  xlab("Year") +
  ylab("Total burn area")
ggsave("fig/ppc-tot-burn-area.pdf", width = 20, height = 10)



# Visualize posterior predictions ---------------------------------
st_y <- m_fit %>%
  rstan::extract() %>%
  `[[`("log_lambda_new") %>%
  reshape2::melt(varnames = c("iter", "reg", "year")) %>%
  group_by(year, reg) %>%
  summarize(median = median(value),
            lo = quantile(value, 0.025),
            hi = quantile(value, 0.975)) %>%
  mutate(par = "log_lambda") %>%
  full_join(data_summary)

st_y %>%
  ggplot(aes(x = fire_year, y = median, group = us_l3name)) +
  geom_line() +
  xlab("Year") +
  ylab("log(Expected fire density)")


st_z <- m_fit %>%
  rstan::extract() %>%
  `[[`("mu_z_new") %>%
  reshape2::melt(varnames = c("iter", "reg", "year")) %>%
  group_by(year, reg) %>%
  summarize(median = median(value),
            lo = quantile(value, 0.025),
            hi = quantile(value, 0.975)) %>%
  mutate(par = "mu_z") %>%
  full_join(data_summary)

st_z %>%
  ggplot(aes(x = fire_year, y = median, group = us_l3name)) +
  geom_line() +
  xlab("Year") +
  ylab("log(Expected fire size)")



jd <- full_join(st_y, st_z) %>%
  dplyr::select(-lo, -hi) %>%
  spread(par, median) %>%
  ungroup %>%
  mutate(year = year + min_year - 1,
        region = regions[reg])

jd %>%
  ggplot(aes(x = log_lambda, y = mu_z, group = reg, color = fire_year)) +
  geom_path() +
  scale_color_viridis() +
  theme_minimal() +
  facet_wrap(~ region) +
  ylab("log(Expected fire size)") +
  xlab("log(Expected fire density)")


scatterD3(jd$log_lambda, jd$mu_z, col_var = jd$region,
          lab = as.character(jd$year),
          xlab = "log(E(number of fires))",
          ylab = "E(log(fire size))")

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
  geom_line(alpha = .8) +
  theme_minimal() +
  scale_y_log10() +
  geom_ribbon(aes(ymin = exp(lo), ymax = exp(hi)), alpha = .03) +
  xlab("Year") +
  ylab("Spatiotemporal adjustment: fire density")

st_y %>%
  full_join(data_summary) %>%
  ungroup() %>%
  ggplot(aes(x = exp(median), y = n_fires / (area))) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  xlab("Predicted fire density") +
  ylab("Observed fire density") +
  scale_x_log10() +
  scale_y_log10()


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
  geom_ribbon(aes(ymin = exp(lo), ymax = exp(hi)), alpha = .03) +
  theme_minimal() +
  xlab("Year") +
  ylab("Spatiotemporal adjustment: burn area")


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
  mutate(us_l3name = regions[reg]) %>%
  left_join(all_ers) %>%
  ggplot(aes(y = reorder(us_l3name, median), x = median, color = over_zero)) +
  geom_point() +
  geom_segment(aes(x = lo, xend = hi, yend = reorder(us_l3name, median))) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  ylab("Ecoregion") +
  xlab("Relationship between fire density and mean burn area") +
  scale_color_brewer(palette = "Set1")
ggsave("fig/n-y-coupling.pdf", width = 8, height = 10)



beta_phi_y %>%
  mutate(us_l3name = regions[reg]) %>%
  left_join(all_ers) %>%
  full_join(er_df) %>%
  ggplot(aes(x = long, y = lat, group = group)) +
  geom_polygon(aes(fill = median), color = NA) +
  coord_equal() +
  labs(x = "Longitude", y = "Latitude") +
  scale_fill_viridis("Relationship between\nfire density and\nmean burn area") +
  theme_map
ggsave("fig/n-y-coupling-map.pdf", width = 12, height = 7)


beta_phi_y %>%
  mutate(us_l3name = regions[reg]) %>%
  left_join(all_ers) %>%
  full_join(er_df) %>%
  ggplot(aes(x = long, y = lat, group = group)) +
  geom_polygon(aes(fill = hi - lo), color = NA) +
  geom_polygon(fill = NA, color = "black", alpha = .1, size = .1) +
  coord_equal() +
  labs(x = "Longitude", y = "Latitude") +
  scale_fill_viridis("Uncertainty in\ndensity-burn area\nrelationsihp",
                     option = "A") +
  theme_map



# Plot spatiotemporal effects  --------------------------------------------
m_fit %>%
  rstan::extract() %>%
  `[[`("phi_y") %>%
  reshape2::melt(varnames = c("iter", "reg", "year")) %>%
  group_by(year, reg) %>%
  summarize(median = median(value)) %>%
  mutate(us_l3name = regions[reg],
         fire_year = year + min_year - 1) %>%
  ungroup %>%
  full_join(er_df) %>%
  ggplot(aes(x = long, y = lat, group = group)) +
  geom_polygon(aes(fill = median), color = NA) +
  coord_equal() +
  labs(x = "Longitude", y = "Latitude") +
  scale_fill_viridis("Expected fire density") +
  theme_map +
  facet_wrap(~ fire_year)
ggsave("fig/n-spatiotemporal.png", width = 20, height = 15)


m_fit %>%
  rstan::extract() %>%
  `[[`("mu_z_new") %>%
  reshape2::melt(varnames = c("iter", "reg", "year")) %>%
  group_by(year, reg) %>%
  summarize(median = median(value)) %>%
  mutate(us_l3name = regions[reg],
         fire_year = year + min_year - 1) %>%
  ungroup %>%
  full_join(er_df) %>%
  ggplot(aes(x = long, y = lat, group = group)) +
  geom_polygon(aes(fill = median), color = NA) +
  coord_equal() +
  labs(x = "Longitude", y = "Latitude") +
  scale_fill_viridis("Expected fire size") +
  theme_map +
  facet_wrap(~ fire_year)
ggsave("fig/y-spatiotemporal.png", width = 20, height = 15)



# Create animated GIF of results ------------------------------------------
years <- unique(d$fire_year) %>% sort
post_counts <- m_fit %>%
  rstan::extract() %>%
  `[[`("phi_y") %>%
  reshape2::melt(varnames = c("iter", "reg", "year")) %>%
  group_by(year, reg) %>%
  summarize(median = median(value)) %>%
  mutate(us_l3name = regions[reg],
         fire_year = year + min_year - 1) %>%
  ungroup %>%
  full_join(er_df)

for (i in seq_along(years)) {
  p <- post_counts %>%
    filter(year == i) %>%
    ggplot(aes(x = long, y = lat, group = group)) +
    geom_polygon(aes(fill = median), color = NA) +
    coord_equal() +
    labs(x = "Longitude", y = "Latitude") +
    scale_fill_viridis("Expected fire density",
                       limits = c(min(post_counts$median),
                                  max(post_counts$median))) +
    theme_map +
    ggtitle(paste(years[i]))
  ggsave(filename = paste0("density-", years[i], ".png"), plot = p)
}





post_size <- m_fit %>%
  rstan::extract() %>%
  `[[`("mu_z_new") %>%
  reshape2::melt(varnames = c("iter", "reg", "year")) %>%
  group_by(year, reg) %>%
  summarize(median = median(value)) %>%
  mutate(us_l3name = regions[reg],
         fire_year = year + min_year - 1) %>%
  ungroup %>%
  full_join(er_df)


for (i in seq_along(years)) {
  p <- post_size %>%
    filter(year == i) %>%
    ggplot(aes(x = long, y = lat, group = group)) +
    geom_polygon(aes(fill = median), color = NA) +
    coord_equal() +
    labs(x = "Longitude", y = "Latitude") +
    scale_fill_viridis("Expected fire size",
                       limits = c(min(post_size$median),
                                  max(post_size$median))) +
    theme_map +
    ggtitle(paste(years[i]))
  ggsave(filename = paste0("size-", years[i], ".png"), plot = p)
}

