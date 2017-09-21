library(tidyverse)
library(rgdal)
library(maptools)
library(foreign)
library(purrr)
library(rasterVis)
library(clusterGeneration)
library(sf)
library(zoo)

source("R/fetch-fire-data.R")

aea_proj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

# Read ecoregion data
ecoregions <- st_read('data/raw/us_eco_l3/us_eco_l3.shp')

# Read fire data ----------------------
mtbs <- st_read('data/raw/mtbs_fod_pts_data/mtbs_fod_pts_20170501.shp') %>%
  filter(!(STATE %in% c("Alaska", "Hawaii", "Puerto Rico")),
         R_ACRES > 1e3 # consistent cutoff for west and east
         ) %>%
  st_transform(st_crs(ecoregions)) %>%
  mutate(ym = as.yearmon(paste(FIRE_YEAR, sprintf("%02d", FIRE_MON),
                               sep = "-")))




# match each ignition to an ecoregion
if (!file.exists("ov.rds")) {
  st_over <- function(x, y) {
    sapply(st_intersects(x,y), function(z) if (length(z)==0) NA_integer_ else z[1])
  }
  ov <- st_over(mtbs, ecoregions)
  write_rds(ov, "ov.rds")
}
ov <- read_rds("ov.rds")

mtbs <- mtbs %>%
  mutate(NA_L3NAME = ecoregions$NA_L3NAME[ov],
         NA_L3NAME = factor(NA_L3NAME, levels = levels(ecoregions$NA_L3NAME))) %>%
  filter(!is.na(NA_L3NAME))

# count the number of fires in each ecoregion in each month
count_df <- mtbs %>%
  tbl_df %>%
  dplyr::select(-geometry) %>%
  group_by(NA_L3NAME, ym) %>%
  summarize(n_fire = n()) %>%
  ungroup %>%
  complete(ym, NA_L3NAME, fill = list(n_fire = 0)) %>%
  arrange(ym)

stopifnot(0 == sum(is.na(count_df$NA_L3NAME)))
stopifnot(sum(count_df$n_fire) == nrow(mtbs))

# load covariate data and link to count data frame
ecoregion_summaries <- read_csv('https://s3-us-west-2.amazonaws.com/earthlab-gridmet/ecoregion_summaries.csv',
                                col_types = cols(
                                  NA_L3NAME = col_character(),
                                  variable = col_character(),
                                  year = col_number(),
                                  month = col_number(),
                                  wmean = col_number())
                                ) %>%
  mutate(year = ifelse(year == 2, 2000, year),
         year = parse_number(year),
         ym = as.yearmon(paste(year,
                 sprintf("%02d", month),
                 sep = "-"))) %>%
  spread(variable, wmean) %>%
  filter(year < 2017)

count_df <- left_join(count_df, ecoregion_summaries)

