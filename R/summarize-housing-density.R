
# Summarizing housing density at the ecoregion level ----------------------
library(raster)
library(snowfall)
library(rgdal)
library(purrr)

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
  dplyr::select(-starts_with("ID")) %>%
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
  ),
  NA_L3NAME = as.character(NA_L3NAME),
  NA_L3NAME = ifelse(NA_L3NAME == 'Chihuahuan Desert',
                     'Chihuahuan Deserts',
                     NA_L3NAME)) %>%
  group_by(NA_L3NAME, year) %>%
  summarize(wmean = weighted.mean(value, Shape_Area)) %>%
  ungroup

# Then interpolate for each month and year from 1984 - 2015
# using a simple linear sequence
impute_density <- function(df) {
  year_seq <- min(df$year):max(df$year)
  predict_seq <- seq(min(df$year),
                     max(df$year),
                     length.out = (length(year_seq) - 1) * 12)
  preds <- approx(x = df$year,
         y = df$wmean,
         xout = predict_seq)
  res <- as_tibble(preds) %>%
    rename(t = x, wmean = y) %>%
    mutate(year = floor(t),
           month = rep(1:12, times = length(year_seq) - 1)) %>%
    filter(year < 2030)
  res$NA_L3NAME <- unique(df$NA_L3NAME)
  res
}

res <- extraction_df %>%
  split(.$NA_L3NAME) %>%
  map(~impute_density(.)) %>%
  bind_rows %>%
  rename(housing_density = wmean)


# evaluate the linear interpolation
extraction_df %>%
  ggplot(aes(year, wmean)) +
  geom_point() +
  geom_point(aes(t, housing_density),
             data = res, color = 'blue',
             size = .1, alpha = .1) +
  facet_wrap(~ NA_L3NAME)


res %>%
  write_csv('data/processed/housing_density.csv')

# system('aws s3 cp data/processed/housing_density.csv s3://earthlab-mjoseph/demo_evt/housing_density.csv --acl public-read')
