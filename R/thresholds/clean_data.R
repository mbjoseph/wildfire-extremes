library(rgdal)
library(dplyr)
library(foreign)
library(purrr)
library(tibble)

ecoregions <- readOGR("data/us_eco_l4", 'us_eco_l4_no_st')
ecoregions <- spTransform(ecoregions, CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'))

d <- read.dbf('data/mtbs/ShortWmtbs.dbf', as.is = TRUE)

d <- d %>%
  filter(!(STATE %in% c('AK', 'PR', 'HI')))
#  filter(STATE == 'ID')

spd <- SpatialPointsDataFrame(coords = d[, c('LONGITUDE', 'LATITUDE')],
                              data = d,
                              proj4string = CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'))

# get ecoregions for each point and merge into data frame
d <- over(spd, ecoregions) %>%
  cbind(d) %>%
  map_if(is.factor, as.character) %>%
  as_tibble()

names(d) <- tolower(names(d))
