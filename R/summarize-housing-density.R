
# Summarizing housing density at the ecoregion level ----------------------
library(raster)
library(snowfall)
library(rgdal)

source('R/helpers.R')

ecoregion_shp <- load_ecoregions()

tifs <- list.files("data/processed",
                   pattern = "*.tif",
                   full.names = TRUE)

extract_one <- function(filename, ecoregion_shp) {
  out_name <- gsub('.tif', '.csv', filename)
  if (!file.exists(out_name)) {
    r <- raster::raster(filename)
    raster::values(r)[raster::values(r) == -999] <- NA
    res <- raster::extract(r, ecoregion_shp,
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

stopifnot(all(lapply(extractions, nrow) == nrow(ecoregion_shp)))

library(tidyverse)

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
  mutate(year = case_when(
    .$variable == 'den00' ~ 2000,
    .$variable == 'den10' ~ 2010,
    .$variable == 'den20' ~ 2020,
    .$variable == 'den30' ~ 2030,
    .$variable == 'den80' ~ 1980,
    .$variable == 'den90' ~ 1990
  )) %>%
  group_by(NA_L3NAME, year) %>%
  summarize(wmean = weighted.mean(value, Shape_Area)) %>%
  ungroup

library(plotly)

p <- extraction_df %>%
  ggplot(aes(year, wmean, group = NA_L3NAME)) +
  geom_line() +
  geom_point() +
  ylab('Average housing density')

ggplotly(p)
