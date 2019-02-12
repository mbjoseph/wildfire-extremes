library(raster)
library(sf)
library(tidyverse)
library(lubridate)
library(patchwork)

# Wallow fire case study  -------------------------------------------------
if (!file.exists('data/raw/mtbs_perimeter_data.zip')) {
  download.file('https://edcintl.cr.usgs.gov/downloads/sciweb1/shared/MTBS_Fire/data/composite_data/burned_area_extent_shapefile/mtbs_perimeter_data.zip' ,
                destfile = 'data/raw/mtbs_perimeter_data.zip')
  dir.create('data/raw/mtbs_perimeters')
  unzip('data/raw/mtbs_perimeter_data.zip', exdir = 'data/raw/mtbs_perimeters/')
}

perimeters <- st_read('data/raw/mtbs_perimeters/mtbs_perims_DD.shp')

wallow_perim <- perimeters %>%
  filter(Fire_Name == 'WALLOW', Year == 2011)

ecoregions <- read_rds('data/processed/ecoregions.rds')

region_color <- 'dodgerblue'


# Evaluate correspondence between local and regional conditions -----------

# load may-july 2011 monthly climate summaries
pr <- raster('data/processed/climate-data/pr/monthly_pr_2011.tif')
rmin <- raster('data/processed/climate-data/rmin/monthly_rmin_2011.tif')
tmmx <- raster('data/processed/climate-data/tmmx/monthly_tmmx_2011.tif')
vs <- raster('data/processed/climate-data/vs/monthly_vs_2011.tif')


# get daily data
read.csv('data/raw/climate-data.csv', stringsAsFactors = FALSE) %>%
  filter(grepl("_2011.nc", url)) %>%
  mutate(basename = basename(url), 
         destfile = file.path('data', 'raw', basename)) %>%
  split(.$destfile) %>%
  map(~download.file(.$url, destfile = .$destfile))

# extract values for the Wallow fire
start_date <- as.Date(paste(2011, "01", "01", sep = "-"))
end_date <- as.Date(paste(2011, "12", "31", sep = "-"))
date_seq <- seq(start_date, end_date, by = "1 day")
is_date_relevant <- date_seq >= as.Date('2011-05-01') & date_seq <= as.Date('2011-07-31')
month_seq <- lubridate::month(date_seq)

daily_pr <- stack('data/raw/pr_2011.nc') %>%
  subset(which(is_date_relevant))
daily_rmin <- stack('data/raw/rmin_2011.nc') %>%
  subset(which(is_date_relevant))
daily_tmmx <- stack('data/raw/tmmx_2011.nc') %>%
  subset(which(is_date_relevant))
daily_vs <- stack('data/raw/vs_2011.nc') %>%
  subset(which(is_date_relevant))

month_seq <- lubridate::month(date_seq)
monthly_pr <- stack('data/raw/pr_2011.nc') %>%
  raster::stackApply(month_seq, fun = mean)
monthly_rmin <- stack('data/raw/rmin_2011.nc') %>%
  raster::stackApply(month_seq, fun = mean)
monthly_tmmx <- stack('data/raw/tmmx_2011.nc') %>%
  raster::stackApply(month_seq, fun = mean)
monthly_vs <- stack('data/raw/vs_2011.nc') %>%
  raster::stackApply(month_seq, fun = mean)

wallow_perim_proj <- st_transform(wallow_perim, projection(daily_pr))

region_proj <- ecoregions %>%
  filter(NA_L3NAME == 'Arizona/New Mexico Mountains') %>%
  st_transform(projection(daily_pr))

local_pr <- raster::extract(daily_pr, as(wallow_perim_proj, 'Spatial'))
local_rmin <- raster::extract(daily_rmin, as(wallow_perim_proj, 'Spatial'))
local_tmmx <- raster::extract(daily_tmmx, as(wallow_perim_proj, 'Spatial'))
local_vs <- raster::extract(daily_vs, as(wallow_perim_proj, 'Spatial'))

region_pr <- raster::extract(monthly_pr, as(region_proj, 'Spatial'), fun = mean) %>%
  apply(2, FUN = function(x) weighted.mean(x, w = region_proj$Shape_Area))
region_rmin <- raster::extract(monthly_rmin, as(region_proj, 'Spatial'), fun = mean) %>%
  apply(2, FUN = function(x) weighted.mean(x, w = region_proj$Shape_Area))
region_tmmx <- raster::extract(monthly_tmmx, as(region_proj, 'Spatial'), fun = mean) %>%
  apply(2, FUN = function(x) weighted.mean(x, w = region_proj$Shape_Area))
region_vs <- raster::extract(monthly_vs, as(region_proj, 'Spatial'), fun = mean) %>%
  apply(2, FUN = function(x) weighted.mean(x, w = region_proj$Shape_Area))


name_df <- tibble(var = c('pr', 'rmin', 'tmmx', 'vs'), 
                  name = c('Precipitation (mm)',
                           'Humidity (%)', 
                           'Maximum air temperature (C)', 
                           'Wind speed (m/s)'))
to_df <- function(mat) {
  colnames(mat) <- as.character(date_seq[is_date_relevant])
  df <- as_tibble(mat) %>%
    mutate(pixel = 1:n()) %>%
    gather(date, value, -pixel)
  return(df)
}

daily_df <- list('pr' = local_pr[[1]], 
                 'rmin' = local_rmin[[1]], 
                 'tmmx' = local_tmmx[[1]] - 273.15, 
                 'vs' = local_vs[[1]]) %>%
  lapply(to_df) %>%
  bind_rows(.id = 'var') %>%
  mutate(date = as.Date(date), 
         month = month(date)) %>%
  left_join(name_df)

region_df <- tibble(pr = region_pr, 
                    rmin = region_rmin, 
                    tmmx = region_tmmx - 273.15, 
                    vs = region_vs, 
                    month = 1:12) %>%
  right_join(distinct(daily_df, date, month)) %>%
  gather(var, value, -month, -date) %>%
  left_join(name_df) %>%
  mutate(label = 'Monthly regional mean')

start_stop_df <- tibble(date = as.Date(c('2011-05-29', '2011-07-08')))

annotation_df <- tibble(date = as.Date(c('2011-05-25', '2011-07-13'))) %>%
  mutate(name = 'Humidity (%)', 
         value = c(50, 55), 
         label = c('Ignited\nMay 29', 'Contained'))

p2 <- daily_df %>%
  mutate(label = 'Daily local conditions') %>%
  ggplot(aes(as.Date(date), value)) + 
  geom_line(data = region_df, 
            aes(group = month, color = label), 
            size = 1) + 
  geom_point(aes(color = label), alpha = .02, size = .6) + 
  facet_wrap(~name, scales = 'free_y', ncol = 1) + 
  theme_minimal() + 
  xlab('') + 
  ylab('') + 
  theme(panel.grid.minor = element_blank(), legend.position = 'top') + 
  geom_vline(aes(xintercept = date), data = start_stop_df, linetype = 'dashed') + 
  scale_x_date(date_labels = "%Y-%m-%d") + 
  geom_text(data = annotation_df, aes(label = label), size = 2.2, 
            hjust = 'right') + 
  scale_color_manual(values = c('red', region_color), 
                     guide = guide_legend(override.aes = list(
                       linetype = c("blank", "solid"),
                       shape = c(19, NA), 
                       alpha = c(1, 1))), 
                     "")
p2
ggsave('fig/wallow-local-conditions.png', width = 5, height = 5)
