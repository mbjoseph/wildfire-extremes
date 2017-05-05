library(raster)

source("R/fetch-cmip6-landuse.R")
# Extracting percent urban land cover from L3 ecoregions ------------------
r <- brick(dest, varname = "urban")
r_years <- 850:2015

source("R/01-clean_data.R")
ecoregion_shp <- readOGR(ecoregion_prefix, "us_eco_l3")
#
# # subset raster to years in data
# r <- subset(r, which(r_years > lower_year & r_years < upper_year))
# r <- projectRaster(r, crs = CRS(proj4string(ecoregion_shp)))
#
# urban_r <- raster::extract(r, ecoregion_shp, df = TRUE)
#
# joined_urban <- ecoregion_shp %>%
#   as.data.frame %>%
#   tbl_df %>%
#   mutate(ID = 1:n()) %>%
#   right_join(urban_r) %>%
#   dplyr::select(US_L3NAME, starts_with("year")) %>%
#   gather(year, urban, -US_L3NAME) %>%
#   mutate(fire_year = substring(year, first = 6, last = 10),
#          fire_year = parse_number(fire_year)) %>%
#   group_by(US_L3NAME, fire_year) %>%
#   summarize(p_urban = mean(urban, na.rm = TRUE)) %>%
#   ungroup
#
# names(joined_urban) <- tolower(names(joined_urban))
#
# write_csv("data/processed/cmip_urban.csv", x = joined_urban)


r <- brick(dest, varname = "secmb")
r_years <- 850:2015

source("R/01-clean_data.R")
ecoregion_shp <- readOGR(ecoregion_prefix, "us_eco_l3")

# subset raster to years in data
r <- subset(r, which(r_years > lower_year & r_years < upper_year))
r <- projectRaster(r, crs = CRS(proj4string(ecoregion_shp)))

secb_r <- raster::extract(r, ecoregion_shp, df = TRUE)

joined_secb <- ecoregion_shp %>%
  as.data.frame %>%
  tbl_df %>%
  mutate(ID = 1:n()) %>%
  right_join(secb_r) %>%
  dplyr::select(US_L3NAME, starts_with("X")) %>%
  gather(year, secb, -US_L3NAME) %>%
  mutate(fire_year = r_years[as.numeric(substring(.$year, first = 2, last = 6))]) %>%
  group_by(US_L3NAME, fire_year) %>%
  summarize(mean_secb = mean(secb, na.rm = TRUE)) %>%
  ungroup

names(joined_secb) <- tolower(names(joined_secb))

write_csv("data/processed/cmip_secb.csv", x = joined_secb)
