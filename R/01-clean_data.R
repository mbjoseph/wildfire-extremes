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
mtbs <- readOGR(mtbs_prefix, 'mtbs_fod_pts_20160401')

explore_short <- FALSE
if (explore_short) {
  # Exploration to explore whether to merge mtbs & short data ----------
  mtbs_df <- mtbs %>%
    as.data.frame() %>%
    tbl_df
  names(mtbs_df) <- tolower(names(mtbs_df))
  mtbs_df <- mtbs_df %>%
    rename(mtbs_id = fire_id,
           mtbs_fire_name = firename,
           latitude = lat,
           longitude = long,
           fire_size = r_acres)

  short <- dbConnect(RSQLite::SQLite(), short_sqlite) %>%
    dbGetQuery('select * from fires') %>%
    dplyr::select(-shape) %>%
    tbl_df

  # select relevant columns from each
  mtbs_df <- mtbs_df %>%
    select(mtbs_id, mtbs_fire_name, fire_year, fire_size, longitude, latitude) %>%
    mutate(source = "mtbs")

  short <- short %>%
    mutate(source = "short") %>%
    select(objectid, fire_name, mtbs_id, mtbs_fire_name, fire_size, fire_year,
           source, longitude, latitude)

  full_join(mtbs_df, short) %>%
    filter(fire_year >= min(short$fire_year),
           fire_size > 1000,
           !duplicated(mtbs_id)) %>%
    ggplot(aes(fire_size, color = source)) +
    geom_density() +
    scale_x_log10()
  # Looks like short hardly adds any large fires - they're mostly on the
  # smaller side. If we care about extremes and longer term temporal patterns,
  # it may be worth using MTBS only.
}

## load EPA level 3 ecoregion data ------------------------------------------
ecoregions <- readOGR(ecoregion_prefix, "us_eco_l3")

# Simplify the polygons for faster plotting & reproject
ecoregions <- ecoregions %>%
  rgeos::gSimplify(tol = 1e3, topologyPreserve = TRUE) %>%
  SpatialPolygonsDataFrame(ecoregions@data)

ecoregions$id <- row.names(ecoregions)
er_df <- fortify(ecoregions, region = 'id')
er_df <- left_join(er_df, ecoregions@data, by = 'id') %>%
  as_tibble()
names(er_df) <- tolower(names(er_df))

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

# subsetting for testing
d <- d %>%
  filter(fire_year > 1989, fire_year < 2000)

lower_year <- min(d$fire_year) - 1
upper_year <- max(d$fire_year) + 1

# clean up
rm(list = c("mtbs", "explore_short"))
