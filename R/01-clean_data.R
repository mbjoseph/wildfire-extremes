library(raster)
library(tidyverse)
library(rgdal)
library(maptools)
library(ggmap)
library(foreign)
library(purrr)
library(RSQLite)
library(rasterVis)
library(clusterGeneration)

source("R/00-fetch-data.R")

# Read fire data ----------------------
mtbs <- readOGR(mtbs_prefix, 'mtbs_fod_pts_20170501') %>%
  subset(!(STATE %in% c("Alaska", "Hawaii", "Puero Rico")))

# Load MACA climate data
source("R/process-maca.R")

aea_proj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "

coarse_rp <- projectRaster(coarse_r, crs = CRS(aea_proj))
mtbs <- spTransform(mtbs, CRS(aea_proj))

projection(coarse_rp)
projection(mtbs)

plot(coarse_rp[[1]])
plot(mtbs, pch = 19, col = 2, add = TRUE, cex = .5)
title("Coarse grid with underlying points")

n_year <- length(unique(mtbs$FIRE_YEAR))
fire_years <- 1984:2015
k <- 1
count_l <- list()
for (i in 1:n_year) {
  count_l[[k]] <- mtbs %>%
    subset(FIRE_YEAR == fire_years[i]) %>%
    rasterize(coarse_rp[[1]], fun = "count", background = 0) %>%
    subset(1)
  k <- k + 1
}

count_r <- stack(count_l)
names(count_r) <- paste0("n_fires_in_", fire_years)

plot(count_r)

count_r <- mask(count_r, coarse_rp[[1]])

plot(count_r)

gplot(count_r) +
  theme_classic() +
  geom_tile(aes(fill = log(value + 1))) +
  facet_wrap(~ variable) +
  scale_fill_viridis("log(# Fires + 1)") +
  coord_equal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggsave(filename = "fig/us-fire-counts.pdf", width = 8, height = 6)


plot(coarse_rp)

count_df <- as.data.frame(count_r) %>%
  tbl_df %>%
  mutate(pixel_idx = 1:n()) %>%
  gather(year, n_fires, -pixel_idx) %>%
  mutate(year = gsub("n_fires_in_", "", year) %>% parse_number)

count_df %>%
  ggplot(aes(year, n_fires, group = pixel_idx)) +
  geom_line()

count_df %>%
  ggplot(aes(n_fires)) +
  geom_histogram()

# verify that the mask is constant over time (this raises error if not)
count_df %>%
  group_by(pixel_idx) %>%
  summarize(p_na = mean(is.na(n_fires)))

# verify that all fires were counted
raster_counts <- count_df %>%
  group_by(year) %>%
  summarize(n_fires = sum(n_fires, na.rm = TRUE))

mtbs_counts <- as.data.frame(mtbs) %>%
  group_by(FIRE_YEAR) %>%
  summarize(n_fires_mtbs = n()) %>%
  rename(year = FIRE_YEAR)

full_join(raster_counts, mtbs_counts) %>%
  ggplot(aes(n_fires_mtbs, n_fires)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0)

## Find the pixel index for each fire
mtbs_covs <- raster::extract(coarse_rp, mtbs, cellnumbers = TRUE)  %>%
  tbl_df %>%
  mutate(mtbs_idx = 1:n()) %>%
  rename(pixel_idx = cells) %>%
  gather(variable, value, -pixel_idx, -mtbs_idx) %>%
  separate(variable, c("variable", "year"), sep = "_") %>%
  spread(variable, value) %>%
  arrange(year, mtbs_idx)


d <- mtbs_covs %>%
  distinct(mtbs_idx, pixel_idx) %>%
  left_join(tbl_df(mtbs) %>% mutate(mtbs_idx = 1:n())) %>%
  map_if(is.factor, as.character) %>%
  as_tibble() %>%
  filter(R_ACRES > 1000) %>%
  mutate(fire_size = R_ACRES)

names(d) <- tolower(names(d))

lower_year <- min(d$fire_year) - 1
upper_year <- max(d$fire_year) + 1

# clean up
rm(list = c("mtbs"))
