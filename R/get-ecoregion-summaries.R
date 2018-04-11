library(tidyverse)
library(raster)
library(snowfall)
library(rgdal)
source('R/helpers.R')


# Extracting monthly climate summaries for ecoregions ---------------------
ecoregion_shp <- rgdal::readOGR(dsn = 'data/raw/us_eco_l3', layer = 'us_eco_l3') %>%
  spTransform(CRSobj = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))

tifs <- list.files("data/processed",
                   pattern = ".tif",
                   recursive = TRUE,
                   full.names = TRUE)

# generate indices for mapping polygons to raster cells
r <- raster::brick(tifs[1])
shp_raster_idx <- cellFromPolygon(r, ecoregion_shp)
names(shp_raster_idx) <- 1:length(shp_raster_idx)

# discard polygons with all NA values in raster
nonzero_list_elements <- shp_raster_idx[!sapply(shp_raster_idx, FUN = is.null)]

# define an efficient extraction function to get mean values by polygon
fast_extract <- function(rasterfile, index_list) {
  r <- raster::brick(rasterfile)
  out <- lapply(index_list, FUN = function(x) {
    raster::extract(r, x) %>%
      colMeans(na.rm = TRUE)
  })
  lapply(out, FUN = function(x) {
    as.data.frame(x) %>%
      rownames_to_column
    }) %>%
    bind_rows(.id = 'ID') %>%
    spread(rowname, x) %>%
    as_tibble
}


# do the extraction
system.time({
  sfInit(parallel = TRUE, cpus = parallel::detectCores())
  sfLibrary(tidyverse)
  extractions <- sfLapply(tifs,
                          fun = fast_extract,
                          index_list = nonzero_list_elements)
  sfStop()
})





# Process extracted values into a usable data frame -----------------------
extraction_df <- extractions %>%
  bind_cols %>%
  mutate(index = ID) %>%
  dplyr::select(-starts_with("ID")) %>%
  rename(ID = index) %>%
  ## TODO: match NA_L3NAME and Shape_Area using the ID column
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
