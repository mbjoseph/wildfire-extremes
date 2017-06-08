library(raster)
library(tidyverse)
library(rgdal)
library(maptools)
library(ggmap)
library(foreign)
library(purrr)
library(RSQLite)


source("R/00-fetch-data.R")

# Read fire data ----------------------
mtbs <- readOGR(mtbs_prefix, 'mtbs_fod_pts_20170501')

# Load MACA climate data
source("R/process-maca.R")

proj_grid <- spTransform(poly_grid, CRS(projection(mtbs)))

plot(proj_grid, col = "grey")
plot(mtbs, pch = 19, col = 2, add = TRUE, cex = .5)

n_year <- length(unique(mtbs$FIRE_YEAR))
fire_years <- 1984:2015
k <- 1
count_l <- list()
for (i in 1:n_year) {
  count_l[[k]] <- mtbs %>%
    subset(FIRE_YEAR == fire_years[i]) %>%
    rasterize(coarse_r[[1]], fun = "count", background = 0) %>%
    subset(1)
  k <- k + 1
}

count_r <- stack(count_l)
names(count_r) <- paste0("n_fires_in_", fire_years)

plot(count_r)

count_r <- mask(count_r, coarse_r[[1]])

plot(count_r)

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
  geom_histogram() +
  facet_wrap(~ year) +
  scale_x_log10()

# verify that the mask is constant over time
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


x <- sp::over(mtbs, proj_grid)

testr <- rasterize(mtbs, coarse_r, fun = "count")

maca <- raster::extract(annual_rstack[[1]], proj_er, fun = mean, df = TRUE,
                        na.rm = TRUE, small = FALSE)
mean(is.na(maca$huss_1984))

maca <- raster::extract(annual_rstack, proj_er, fun = mean, df = TRUE)



# overlay mtbs fire data onto ecoregion shapefile
d <- mtbs %>%
  spTransform(CRS(proj4string(ecoregions))) %>%
  over(ecoregions) %>%
  bind_cols(as.data.frame(mtbs)) %>%
  map_if(is.factor, as.character) %>%
  as_tibble() %>%
  filter(!is.na(US_L3NAME),
         R_ACRES > 1000) %>%
  mutate(fire_size = R_ACRES)

names(d) <- tolower(names(d))

# # subsetting for testing
# d <- d %>%
#   filter(fire_year > 1989, fire_year < 2000)

lower_year <- min(d$fire_year) - 1
upper_year <- max(d$fire_year) + 1

# clean up
rm(list = c("mtbs", "explore_short"))
