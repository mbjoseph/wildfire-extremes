library(sf)
library(tidyverse)
library(sp)
library(elevatr)
library(raster)

ecoregions <- st_read('data/raw/us_eco_l3/us_eco_l3.shp') %>%
  as('Spatial') %>%
  spTransform(CRS('+init=epsg:4326'))

elevation <- get_elev_raster(ecoregions, z = 8, src = 'aws')

plot(elevation)
plot(ecoregions, add = TRUE)

tri <- terrain(elevation, opt = c('TRI'))
plot(tri)

ecoregion_tri <- raster::extract(tri, ecoregions, fun = mean, sp = TRUE)

# now take a weighted average, weighting each subpolygon by area
tri_df <- ecoregion_tri %>%
  as.data.frame %>%
  tbl_df %>%
  mutate(NA_L3NAME = as.character(NA_L3NAME),
         NA_L3NAME = ifelse(NA_L3NAME == 'Chihuahuan Desert',
                            'Chihuahuan Deserts',
                            NA_L3NAME)) %>%
  group_by(NA_L3NAME) %>%
  summarize(tri = weighted.mean(tri, Shape_Area)) %>%
  arrange(-tri)

# write to csv and then push to Amazon S3
# with aws s3 cp from command line
tri_df %>%
  write_csv('ecoregion_tri.csv')

system('aws s3 cp ecoregion_tri.csv s3://earthlab-gridmet/ecoregion_tri.csv --acl public-read')
