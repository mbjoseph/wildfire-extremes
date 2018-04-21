library(tidyverse)
library(raster)
library(snowfall)
library(rgdal)
library(assertthat)
source('R/helpers.R')


# Extracting monthly climate summaries for ecoregions ---------------------
ecoregion_shp <- rgdal::readOGR(dsn = 'data/raw/us_eco_l3', layer = 'us_eco_l3') %>%
  spTransform(CRSobj = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))

ecoregion_shp$NA_L3NAME <- as.character(ecoregion_shp$NA_L3NAME)

ecoregion_shp$NA_L3NAME <- ifelse(ecoregion_shp$NA_L3NAME == 'Chihuahuan Desert',
                                  'Chihuahuan Deserts',
                                  ecoregion_shp$NA_L3NAME)

tifs <- list.files("data/processed",
                   pattern = ".tif",
                   recursive = TRUE,
                   full.names = TRUE)

# remove any housing density geotiffs that matched the file listing
tifs <- tifs[!grepl('den[0-9]{2}\\.tif', tifs)]

# Generate indices from polygons for raster extraction --------------------
r <- raster::brick(tifs[1])
shp_raster_idx <- cellFromPolygon(r, ecoregion_shp)
names(shp_raster_idx) <- ecoregion_shp$NA_L3NAME

# this list of indices has one element per polygon, but we want one per region
ecoregion_raster_idx <- vector(mode = 'list',
                               length = length(unique(ecoregion_shp$NA_L3NAME)))
ecoregion_names <- sort(unique(ecoregion_shp$NA_L3NAME))
names(ecoregion_raster_idx) <- ecoregion_names
for (i in seq_along(ecoregion_names)) {
  list_elements <- names(shp_raster_idx) == ecoregion_names[i]
  assert_that(any(list_elements))
  ecoregion_raster_idx[[i]] <- shp_raster_idx[list_elements] %>%
    unlist
}

# verify that no cells are duplicated
assert_that(ecoregion_raster_idx %>%
              sapply(FUN = function(x) any(duplicated(x))) %>%
              sum == 0)

# verify that all ecoregions have some cells
assert_that(ecoregion_raster_idx %>%
              sapply(FUN = function(x) length(x)) %>%
              min > 0)

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
    bind_rows(.id = 'NA_L3NAME') %>%
    spread(rowname, x) %>%
    as_tibble
}




# Extract climate data ---------------------------------------
system.time({
  sfInit(parallel = TRUE, cpus = parallel::detectCores())
  sfLibrary(tidyverse)
  extractions <- sfLapply(tifs,
                          fun = fast_extract,
                          index_list = ecoregion_raster_idx)
  sfStop()
})



# Process extracted values into a usable data frame -----------------------
ecoregion_summaries <- extractions %>%
  bind_cols %>%
  gather(variable, value, -NA_L3NAME) %>%
  filter(!grepl(pattern = 'NA_L3NAME', x = variable)) %>%
  separate(variable,
           into = c("interval", "variable", "timestep"),
           sep = "_") %>%
  separate(timestep, into = c("year", "month"), sep = "\\.") %>%
  dplyr::select(-interval) %>%
  mutate(year = parse_number(year),
         month = parse_number(month)) %>%
  arrange(year, month, variable, NA_L3NAME)

destfile <- "data/processed/ecoregion_summaries.csv"
write_csv(ecoregion_summaries, destfile)

print(paste('Ecoregion climate summaries written to', destfile))
