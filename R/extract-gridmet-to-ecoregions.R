
# Compute ecoregion means of input layers ---------------------------------
library(tidyverse)
library(sf)
library(raster)
library(sp)
library(gdal)
library(geosphere)
library(viridis)

layers <- list.files("data/processed", pattern = "monthly",
                     full.names = TRUE)
r <- stack(layers)

ecoregions <- st_read('data/raw/us_eco_l3/us_eco_l3.shp') %>%
  st_transform(projection(r)) %>%
  as("Spatial")

vals <- raster::extract(r,
                        ecoregions,
                        na.rm = TRUE,
                        fun = mean,
                        df = TRUE)


vals$area_sq_km <- areaPolygon(ecoregions) / 1000

ecoregion_df <- as.data.frame(ecoregions)

vals$US_L3NAME <- ecoregion_df$US_L3NAME

er_vals <- vals %>%
  tbl_df %>%
  gather(var, value, -ID, -area_sq_km, -US_L3NAME) %>%
  filter(!is.na(value)) %>%
  group_by(US_L3NAME, var) %>%
  summarize(mean_val = weighted.mean(value, area_sq_km)) %>%
  arrange(var) %>%
  separate(var, into = c("timestep", "var", "year"), sep = "_") %>%
  separate(year, into = c("year", "month"), sep = "\\.") %>%
  mutate(year = parse_number(year),
         month = parse_number(month)) %>%
  ungroup

er_avgs <- er_vals %>%
  group_by(month) %>%
  summarize(avg = mean(mean_val))

er_vals %>%
  ggplot(aes(month, mean_val, color = year, group = year)) +
  geom_line(data = er_avgs, aes(x = month, y = avg),
            inherit.aes = FALSE) +
  geom_line(alpha = .5) +
  facet_wrap(~US_L3NAME) +
  scale_color_viridis()

write_csv(er_vals, 'data/processed/extracted_pet.csv')
