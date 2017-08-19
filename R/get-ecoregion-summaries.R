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
  raster::extract(raster::stack(filename), ecoregion_shp,
                  na.rm = TRUE, fun = mean, df = TRUE)
}

sfInit(parallel = TRUE, cpus = 40)
sfExport(list = c("ecoregion_shp"))

extractions <- sfLapply(as.list(tifs),
                        fun = extract_one,
                        ecoregion_shp = ecoregion_shp)
sfStop()

# convert to a data frame
library(tidyverse)
ecoregion_summaries <- extractions %>%
  bind_cols %>%
  as_tibble %>%
  mutate(index = ID) %>%
  select(-starts_with("ID")) %>%
  rename(ID = index) %>%
  mutate(NA_L3NAME = data.frame(ecoregion_shp)$NA_L3NAME,
         Shape_Area = data.frame(ecoregion_shp)$Shape_Area) %>%
  gather(variable, value, -NA_L3NAME, -Shape_Area, -ID) %>%
  filter(!is.na(value)) %>%
  separate(variable,
           into = c("interval", "variable", "timestep"),
           sep = "_") %>%
  separate(timestep, into = c("year", "month"), sep = "\\.") %>%
  group_by(NA_L3NAME, variable, year, month) %>%
  summarize(wmean = weighted.mean(value, Shape_Area)) %>%
  ungroup %>%
  mutate(year = parse_number(year),
         month = parse_number(month)) %>%
  arrange(year, month, variable, NA_L3NAME)

write_csv(ecoregion_summaries, "data/processed/ecoregion_summaries.csv")
