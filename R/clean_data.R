library(tidyverse)
library(rgdal)
library(maptools)
library(ggmap)
library(foreign)
library(purrr)
library(raster)

# load ecoregion data
level <- 3
ecoregions <- readOGR(paste0("data/raw/us_eco_l", level),
                      paste0("us_eco_l", level))

# Simplify the polygons for faster plotting & reproject ----------------------
wgs84 <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

ecoregions <- ecoregions %>%
  rgeos::gSimplify(tol = 1e3, topologyPreserve = TRUE) %>%
  SpatialPolygonsDataFrame(ecoregions@data) %>%
  spTransform(CRS(wgs84))

# Create a data frame from ecoregion data ---------------------------------
ecoregions$id <- row.names(ecoregions)
er_df <- fortify(ecoregions, region = 'id')
er_df <- left_join(er_df, ecoregions@data, by = 'id') %>%
  as_tibble()
names(er_df) <- tolower(names(er_df))

# Read short and mtbs fire data & merge with ecoregions ----------------------
mtbs <- readOGR('data/raw/mtbs', 'ShortWmtbs') %>%
  subset(!(STATE %in% c("AK", "PR", "HI"))) %>%
  spTransform(CRS(wgs84))

d <- mtbs %>%
  over(ecoregions) %>%
  bind_cols(as.data.frame(mtbs)) %>%
  map_if(is.factor, as.character) %>%
  as_tibble() %>%
  filter(!is.na(US_L3NAME))

names(d) <- tolower(names(d))

## Fetch pdsi at fire locations
source("R/read_pdsi.R")

pdsi <- raster::extract(pdsi,
                                 y = ecoregions, fun = mean, df = TRUE, sp = TRUE) %>%
  tbl_df() %>%
  dplyr::select(US_L3NAME, Shape_Area, starts_with("X")) %>%
  gather(key = year, value = mean_pdsi, -US_L3NAME, -Shape_Area) %>%
  mutate(year = gsub("X", "", x = year),
         year = as.numeric(year)) %>%
  group_by(US_L3NAME, year) %>%
  # next line computes a weighted average across all polygons in an ecoregion
  summarize(same_year_pdsi = sum(Shape_Area * mean_pdsi, na.rm = TRUE) /
              sum(Shape_Area, na.rm = TRUE)) %>%
  mutate(fire_year = year) %>%
  dplyr::select(-year) %>%
  ungroup()

names(pdsi) <- tolower(names(pdsi))


## Fetch mean annual precip at fire locations --------------------------

# create a raster brick of mean precip
precip_d <- read_csv("data/processed/cleaned-precip.csv")
precip_raster <- list()
years <- 1991:max(d$fire_year)
for (i in seq_along(years)) {
  spg <- precip_d %>%
    filter(year == years[i]) %>%
    dplyr::select(longitude, latitude, mean_precip)
  coordinates(spg) <- ~longitude + latitude
  gridded(spg) <- TRUE
  precip_raster[[i]] <- raster(spg)
}
names(precip_raster) <- years
precip_raster <- brick(precip_raster)
projection(precip_raster) <- CRS("+init=epsg:4326")
precip_raster <- projectRaster(precip_raster, crs = CRS(wgs84))

# compute the mean of the mean annual precip for each ecoregion X year
extracted_precip <- raster::extract(precip_raster,
                y = ecoregions, fun = mean, df = TRUE, sp = TRUE) %>%
  tbl_df() %>%
  dplyr::select(US_L3NAME, Shape_Area, starts_with("X")) %>%
  gather(key = year, value = mean_precip, -US_L3NAME, -Shape_Area) %>%
  mutate(year = gsub("X", "", x = year),
         year = as.numeric(year)) %>%
  group_by(US_L3NAME, year) %>%
  # next line computes a weighted average across all polygons in an ecoregion
  summarize(same_year_precip = sum(Shape_Area * mean_precip, na.rm = TRUE) /
              sum(Shape_Area, na.rm = TRUE)) %>%
  mutate(fire_year = year) %>%
  dplyr::select(-year) %>%
  ungroup()

names(extracted_precip) <- tolower(names(extracted_precip))

precip <- extracted_precip %>%
  mutate(fire_year = fire_year + 1,
         last_year_precip = same_year_precip) %>%
  dplyr::select(-same_year_precip) %>%
  full_join(extracted_precip) %>%
  filter(fire_year > 1991, fire_year < 2014)

rm(list = c("level", "wgs84", "extracted_precip", "mtbs"))
