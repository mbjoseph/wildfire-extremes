library(raster)
library(snowfall)
library(rgdal)
source('R/helpers.R')

ecoregion_shp <- load_ecoregions()

# simplify the shapefile a bit - it's much finer than the rasters
system('ogr2ogr -progress -simplify 40 data/processed/simple-ecoregions.shp data/raw/us_eco_l3/us_eco_l3.shp')

ecoregion_shp <- rgdal::readOGR(dsn = "data/processed/",
                                layer = "simple-ecoregions")
ecoregion_shp <- spTransform(ecoregion_shp,
                       CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))

tifs <- list.files("data/processed",
                   pattern = ".tif",
                   recursive = TRUE,
                   full.names = TRUE)

extract_one <- function(filename, ecoregion_shp) {
  out_name <- gsub('.tif', '.csv', filename)
  if (!file.exists(out_name)) {
    res <- raster::extract(raster::stack(filename), ecoregion_shp,
                    na.rm = TRUE, fun = mean, df = TRUE)
    write.csv(res, file = out_name)
  } else {
    res <- read.csv(out_name)
  }
  res
}

sfInit(parallel = TRUE, cpus = parallel::detectCores())
sfExport(list = c("ecoregion_shp"))

extractions <- sfLapply(as.list(tifs),
                        fun = extract_one,
                        ecoregion_shp = ecoregion_shp)
sfStop()

# ensure that they all have the same length
stopifnot(all(lapply(extractions, nrow) == nrow(ecoregion_shp)))

library(tidyverse)

# push to S3
write_rds(extractions, 'extractions.rds')
system('aws s3 cp extractions.rds s3://earthlab-mjoseph/demo_evt/extractions.rds')

# convert to a data frame
extraction_df <- extractions %>%
  bind_cols %>%
  as_tibble %>%
  mutate(index = ID) %>%
  select(-starts_with("ID")) %>%
  rename(ID = index) %>%
  mutate(NA_L3NAME = data.frame(ecoregion_shp)$NA_L3NAME,
         Shape_Area = data.frame(ecoregion_shp)$Shape_Area) %>%
  dplyr::select(-starts_with('X')) %>%
  gather(variable, value, -NA_L3NAME, -Shape_Area, -ID) %>%
  filter(!is.na(value)) %>%
  mutate(NA_L3NAME = as.character(NA_L3NAME),
         NA_L3NAME = ifelse(NA_L3NAME == 'Chihuahuan Desert',
                            'Chihuahuan Deserts',
                            NA_L3NAME))

ecoregion_summaries <- extraction_df %>%
  separate(variable,
           into = c("interval", "variable", "timestep"),
           sep = "_") %>%
  separate(timestep, into = c("year", "month"), sep = "\\.") %>%
  select(-interval) %>%
  mutate(NA_L3NAME = as.character(NA_L3NAME),
         NA_L3NAME = ifelse(NA_L3NAME == 'Chihuahuan Desert',
                            'Chihuahuan Deserts',
                            NA_L3NAME)) %>%
  group_by(NA_L3NAME, variable, year, month) %>%
  summarize(wmean = weighted.mean(value, Shape_Area)) %>%
  ungroup %>%
  mutate(year = parse_number(year),
         month = parse_number(month)) %>%
  arrange(year, month, variable, NA_L3NAME)

write_csv(ecoregion_summaries, "ecoregion_summaries.csv")
