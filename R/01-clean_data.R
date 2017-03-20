library(tidyverse)
library(rgdal)
library(maptools)
library(ggmap)
library(foreign)
library(purrr)
library(raster)

# load EPA level 3 ecoregion data
# accessible from:
# https://www.epa.gov/eco-research/level-iii-and-iv-ecoregions-continental-united-states
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
# Data from:
# http://www.mtbs.gov/nationalregional/pointdata.html
mtbs <- readOGR('data/raw/mtbs_fod_pts_data/', 'mtbs_fod_pts_20160401') %>%
  subset(!(STATE %in% c("AK", "PR", "HI"))) %>%
  spTransform(CRS(wgs84))

d <- mtbs %>%
  over(ecoregions) %>%
  bind_cols(as.data.frame(mtbs)) %>%
  map_if(is.factor, as.character) %>%
  as_tibble() %>%
  filter(!is.na(US_L3NAME)) %>%
  mutate(fire_size = R_ACRES)

names(d) <- tolower(names(d))

## Fetch pdsi at fire locations
rm(list = c("level", "wgs84",  "mtbs"))
