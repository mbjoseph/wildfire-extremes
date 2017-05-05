# TODO: put this in a makefile with target cmip_precip.csv

# Extract average precipitation for each ecoregoin  -----------------------
source("R/fetch-cmip5-precip.R")
source("R/01-clean_data.R")
ecoregion_shp <- readOGR(ecoregion_prefix, "us_eco_l3")

prcp_years <- substring(names(rstack), first = 37, last = 40) %>%
  as.numeric


# subset precip data to relevant years
rstack <- subset(rstack, which(prcp_years >= lower_year & prcp_years <= upper_year))
rstack <- projectRaster(rstack, crs = CRS(proj4string(ecoregion_shp)))

prcp_r <- raster::extract(rstack, ecoregion_shp, df = TRUE)

# add a column with the ecoregion name and then group by er and compute mean
joined_prcp <- ecoregion_shp %>%
  as.data.frame %>%
  tbl_df %>%
  mutate(ID = 1:n()) %>%
  right_join(prcp_r) %>%
  dplyr::select(US_L3NAME, starts_with("ncar")) %>%
  gather(year, prcp, -US_L3NAME) %>%
  mutate(fire_year = as.numeric(substring(year, first = 37, last = 40))) %>%
  group_by(US_L3NAME, fire_year) %>%
  summarize(prcp = mean(prcp, na.rm = TRUE)) %>%
  ungroup

names(joined_prcp) <- tolower(names(joined_prcp))

joined_prcp <- joined_prcp %>%
  mutate(prcp_lag1 = prcp[match(fire_year - 1, fire_year)]) %>%
  na.omit

# save as a csv file
write_csv("data/processed/cmip_precip.csv", x = joined_prcp)
