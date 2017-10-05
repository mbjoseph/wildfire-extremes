library(raster)
library(snowfall)
library(rgdal)

if (!dir.exists("data/raw/us_eco_l3/")) {
  download.file("ftp://newftp.epa.gov/EPADataCommons/ORD/Ecoregions/us/us_eco_l3.zip",
                destfile = "data/raw/us_eco_l3.zip")
  unzip("data/raw/us_eco_l3.zip",
        exdir = "data/raw/us_eco_l3/")
}

ecoregion_shp <- readOGR(dsn = "data/raw/us_eco_l3/",
                   layer = "us_eco_l3")

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

# convert to a data frame
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
  filter(!is.na(value))

# pdsi has a different naming scheme, need to process differently
pdsi_df <- extraction_df %>%
  filter(grepl('pdsi', variable)) %>%
  separate(variable,
           into = c('variable', 'year', 'month'),
           sep = "_") %>%
  filter(year > 1982)

non_pdsi_df <- extraction_df %>%
  filter(!grepl('pdsi', variable)) %>%
  separate(variable,
           into = c("interval", "variable", "timestep"),
           sep = "_") %>%
  separate(timestep, into = c("year", "month"), sep = "\\.") %>%
  select(-interval)

ecoregion_summaries <- full_join(pdsi_df, non_pdsi_df) %>%
  group_by(NA_L3NAME, variable, year, month) %>%
  summarize(wmean = weighted.mean(value, Shape_Area)) %>%
  ungroup %>%
  mutate(year = parse_number(year),
         month = parse_number(month)) %>%
  arrange(year, month, variable, NA_L3NAME)

write_csv(ecoregion_summaries, "data/processed/ecoregion_summaries.csv")

