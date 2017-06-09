library(scales)
library(gridExtra)
library(spdep)
library(viridis)
library(tidyverse)
library(lubridate)


source("R/01-clean_data.R")

# Visualize yearly fire size distributions ----------------------------
yd <- d %>%
  group_by(fire_year) %>%
  summarize(mean_size = mean(fire_size),
            max_size = max(fire_size)) %>%
  ungroup %>%
  tidyr::complete(fire_year,
           fill = list(mean_size = NA,
                       max_size = NA)) %>%
  gather(Measure, value, -fire_year)


yd %>%
  ggplot(aes(x = fire_year, y = value)) +
  geom_path() +
  facet_wrap(~ Measure, scales = "free_y") +
  xlab("Time") +
  ylab("Value") +
  theme_minimal()
ggsave("fig/explore/mean-max-ts.pdf")

# Create spatial neighbors ------------------------------------------------
W <- coarse_rp %>%
  rasterToPolygons %>%
  poly2nb %>%
  nb2mat(zero.policy = TRUE, style = 'B')
plot(raster(W))

plot(rasterToPolygons(coarse_rp))
coarse_rp %>%
  rasterToPolygons %>%
  poly2nb %>%
  plot(col = "red",
       coords = coordinates(rasterToPolygons(coarse_rp)))

# Create spatiotemporal covariate data frame
st_cov_chunk <- coarse_rp %>%
  as.data.frame(xy = TRUE) %>%
  tbl_df %>%
  mutate(pixel_idx = 1:n()) %>%
  gather(variable, value, -x, -y, -pixel_idx) %>%
  separate(variable, c("variable", "year"), sep = "_") %>%
  arrange(pixel_idx, year, variable) %>%
  filter(!is.na(value)) %>%
  spread(variable, value) %>%
  arrange(pixel_idx, year)

st_covs <- st_cov_chunk %>%
  mutate(dim = 1) %>%
  full_join(mutate(st_cov_chunk, dim = 2)) %>%
  mutate(chuss = c(scale(huss)),
         cprec = c(scale(pr)),
         ctmax = c(scale(tasmax)),
         ctmin = c(scale(tasmin)),
         cwas = c(scale(was)),
         year = parse_number(year),
         dim = factor(dim))


# Create design matrices for each timestep
X <- array(dim = c(n_year, nrow(st_covs) / n_year, 7)) # 2: bivariate
for (i in 1:n_year) {
  X[i, , ] <- filter(st_covs, year - min(year) + 1 == i) %>%
    model.matrix(~ dim + chuss + cprec + ctmax + ctmin + cwas, data = .)
}

A_t <- bdiag(replicate(2, list(W)))
image(A_t)

# compute multivariate basis vectors
r <- 10 # number of basis vectors

S_X <- array(dim = c(n_year, nrow(st_covs) / n_year, r))
for (i in 1:n_year) {
  G <- (diag(nrow(X[1, , ])) -
          X[i, , ] %*% solve(t(X[i, , ]) %*% X[i, , ]) %*% t(X[i, , ])) %*%
    A_t %*%
    (diag(nrow(X[1, , ])) -
       X[i, , ] %*% solve(t(X[i, , ]) %*% X[i, , ]) %*% t(X[i, , ]))
  eG <- eigen(G)
  basis_vectors <- eG$vectors
  S_X[i, , ] <- basis_vectors[, 1:r]
}

# visualize multivariate basis functions
G_unrolled <- S_X[1, , ]
for (i in 2:n_year) G_unrolled = rbind(G_unrolled, S_X[i, , ])

tbl_df(G_unrolled) %>%
  bind_cols(st_covs) %>%
  select(x, y, year, dim, starts_with("V")) %>%
  gather(basis_dim, value, -x, -y, -year, -dim) %>%
  ggplot(aes(x, y, fill = value)) +
  geom_tile() +
  facet_grid(interaction(basis_dim, dim) ~ year) +
  scale_fill_viridis() +
  coord_equal()




## Add area & num_neighbors to data frame -------------------------------
a_df <- data.frame(id = sapply(slot(ecoregions, "polygons"), slot, "ID"),
                   area = sapply(slot(ecoregions, "polygons"), slot, "area"),
                   US_L3NAME = ecoregions@data$US_L3NAME,
                   stringsAsFactors = FALSE) %>%
  group_by(US_L3NAME) %>%
  summarize(area = sum(area))
a_df$n_neighbors <- rowSums(W)
names(a_df) <- tolower(names(a_df))

er_df <- left_join(er_df, a_df) %>%
  mutate(fire_den = n_fires / area,
         log_fire_density = log(n_fires) - log(area))

# get area of each ecoregion
all_ers <- er_df %>%
  group_by(us_l3name) %>%
  summarize(area = unique(area))



# Visualize number of neighbors --------------------------------------
# theme_map <- theme(axis.line = element_blank(),
#                    axis.text.x = element_blank(),
#                    axis.text.y = element_blank(),
#                    axis.ticks = element_blank(),
#                    axis.title.x = element_blank(),
#                    axis.title.y = element_blank(),
#                    panel.background = element_blank(),
#                    panel.border = element_blank(),
#                    panel.grid.major = element_blank(),
#                    panel.grid.minor = element_blank(),
#                    plot.background = element_blank())

# ggplot(er_df, aes(long, lat, group = group)) +
#   geom_polygon(aes(fill = n_neighbors), color = NA) +
#   coord_equal() +
#   scale_fill_viridis() +
#   theme_map
#
#
# # Visualize the fire density in each ecoregion ----------------------------
# ggplot(er_df, aes(long, lat, group = group)) +
#   geom_polygon(aes(fill = log_fire_density), color = NA) +
#   coord_equal() +
#   labs(x = "Longitude", y = "Latitude") +
#   scale_fill_viridis("log(Fire density)") +
#   theme_map

# get count data (number of fires in each ecoregion X year)
data_summary <- d %>%
  group_by(us_l3name, fire_year, fire_mon) %>%
  summarize(n_fires = n()) %>%
  left_join(a_df) %>%
  full_join(all_ers) %>%
  ungroup() %>%
  dplyr::select(-area, -n_neighbors) %>%
  complete(us_l3name, fire_year, fire_mon,
           fill = list(n_fires = 0)) %>%
  full_join(all_ers) %>%
  mutate(cyear = c(scale(fire_year)),
         year = fire_year + 1 - min(fire_year),
         freg = factor(us_l3name,
                       levels = levels(factor(all_ers$us_l3name))),
         reg = as.numeric(freg),
         num_ym = as.numeric(factor(paste(fire_year,
                                          sprintf("%02d", fire_mon)))))

ymdf <- data_summary %>%
  distinct(fire_year, fire_mon, num_ym) %>%
  mutate(ymd = ymd(paste(fire_year,
                         sprintf("%02d", fire_mon),
                         sprintf("%02d", 1),
                         sep = "-")))


# get fire size data
fire_sizes <- d %>%
  dplyr::select(us_l3name, fire_year, fire_mon, fire_size) %>%
  mutate(cyear = c(scale(fire_year)),
         freg = factor(us_l3name,
                       levels = levels(data_summary$freg)),
         reg = as.numeric(freg),
         year = fire_year + 1 - min(fire_year)) %>%
  arrange(us_l3name, fire_year, fire_mon) %>%
  left_join(ymdf)


data_summary %>%
  ggplot(aes(as.numeric(fire_mon), n_fires / area, color = year)) +
  geom_path(alpha = .5) +
  scale_y_log10() +
  facet_wrap(~ us_l3name) +
  scale_color_viridis() +
  theme_minimal()
ggsave("fig/explore/fire-density-seasonality.pdf")

fire_sizes %>%
  ggplot(aes(as.numeric(fire_mon), fire_size, color = year)) +
  geom_jitter(alpha = .5, width = .3, size = .6) +
  scale_y_log10() +
  facet_wrap(~ us_l3name) +
  scale_color_viridis() +
  theme_minimal()
ggsave("fig/explore/fire-size-seasonality.pdf")

fire_sizes %>%
  ggplot(aes(ymd, fire_size, group = year)) +
  geom_point(alpha = .5) +
  scale_y_log10() +
  facet_wrap(~ us_l3name) +
  xlab("Year") +
  ylab("Fire size")
ggsave("fig/explore/fire-size-ts.pdf")


data_summary %>%
  left_join(ymdf) %>%
  ggplot(aes(ymd, n_fires/area)) +
  geom_line(alpha = .5) +
  scale_y_log10() +
  facet_wrap(~ us_l3name) +
  xlab("Year") +
  ylab("Fire density")
ggsave("fig/explore/fire-density-ts.pdf")
