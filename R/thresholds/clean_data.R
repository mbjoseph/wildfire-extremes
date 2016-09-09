library(rgdal)
library(maptools)
library(dplyr)
library(ggmap)
library(foreign)
library(purrr)
library(tibble)

ecoregions <- readOGR("data/us_eco_l3", "us_eco_l3")

simple_ecoregions = rgeos::gSimplify(ecoregions, tol = 100, topologyPreserve = TRUE)
ecoregions <- SpatialPolygonsDataFrame(simple_ecoregions, ecoregions@data)
#
# # remove channel islands - they have no neighbors, only 2 fires
# ecoregions <- ecoregions[ecoregions$US_L3NAME != 'Southern Channel Islands' &
#                            ecoregions$US_L3NAME != 'Northern Channel Islands', ]
#
# ecoregions@data$US_L3NAME <- droplevels.factor(ecoregions@data$US_L3NAME)

ecoregions <- spTransform(ecoregions, CRS("+proj=laea +lat_0=45.5 +lon_0=-114.125 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"))
ecoregions$id <- row.names(ecoregions)
er_df <- fortify(ecoregions, region = 'id')
er_df <- left_join(er_df, ecoregions@data, by = 'id')
names(er_df) <- tolower(names(er_df))

d <- read.dbf('data/mtbs/ShortWmtbs.dbf', as.is = TRUE)

d <- d %>%
  filter(!(STATE %in% c('AK', 'PR', 'HI')))

spd <- SpatialPointsDataFrame(coords = d[, c('LONGITUDE', 'LATITUDE')],
                              data = d,
                              proj4string = CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'))
spd <- spTransform(spd, CRS(proj4string(ecoregions)))
d <- over(spd, ecoregions) %>%
  cbind(d) %>%
  map_if(is.factor, as.character) %>%
  as_tibble()

names(d) <- tolower(names(d))

d <- subset(d, !is.na(us_l3name))
