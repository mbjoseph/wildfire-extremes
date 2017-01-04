library(rgdal)
library(maptools)
library(dplyr)
library(ggmap)
library(foreign)
library(purrr)
library(tibble)

# load ecoregion data
level <- 3
ecoregions <- readOGR(paste0("data/us_eco_l", level),
                      paste0("us_eco_l", level))

# Simplify the polygons for faster plotting & reproject ----------------------
wgs84 <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
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
d <- read.dbf('data/mtbs/ShortWmtbs.dbf', as.is = TRUE) %>%
  filter(!(STATE %in% c('AK', 'PR', 'HI')))

d <- SpatialPointsDataFrame(coords = d[, c('LONGITUDE', 'LATITUDE')],
                              data = d,
                              proj4string = CRS(wgs84)) %>%
  over(ecoregions) %>%
  cbind(d) %>%
  map_if(is.factor, as.character) %>%
  as_tibble() %>%
  filter(!is.na(US_L3NAME))
names(d) <- tolower(names(d))

rm(list = c("level", "wgs84"))
