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



## Get perimeters for the Wallow Fire from GeoMAC --------------------------------------------
wallow_shp_dir <- 'data/processed/wallow_shp'
dir.create(wallow_shp_dir)
geomac_wallow_prefix <- 'https://rmgsc.cr.usgs.gov/outgoing/GeoMAC/2011_fire_data/Arizona/Wallow/'
wallow_zips <- c(
  'az_wallow_20110530_1130_dd83.zip', 
  'az_wallow_20110531_2332_dd83.zip', 
  'az_wallow_20110602_0013_dd83.zip',
  'az_wallow_20110603_0025_dd83.zip',
  'az_wallow_20110605_2209_dd83.zip',
  'az_wallow_20110606_2304_dd83.zip',
  'az_wallow_20110607_2114_dd83.zip',
  'az_wallow_20110608_2357_dd83.zip',
  'az_wallow_20110609_2013_dd83.zip',
  'az_wallow_20110610_2329_dd83.zip',
  'az_wallow_20110611_2231_dd83.zip',
  'az_wallow_20110612_2316_dd83.zip',
  'az_wallow_20110613_2238_dd83.zip',
  'az_wallow_20110614_2344_dd83.zip',
  'az_wallow_20110615_2220_dd83.zip', 
  'az_wallow_20110616_2330_dd83.zip',
  'az_wallow_20110617_2343_dd83.zip', 
  'az_wallow_20110618_1257_dd83.zip', 
  'az_wallow_20110619_2349_dd83.zip',
  'az_wallow_20110620_2338_dd83.zip',
  'az_wallow_20110621_2345_dd83.zip',
  'az_wallow_20110622_2349_dd83.zip',
  'az_wallow_20110624_0045_dd83.zip',
  'az_wallow_20110626_2323_dd83.zip',
  'az_wallow_20110627_2249_dd83.zip'
)
for (i in seq_along(wallow_zips)) {
  destfile <- file.path(wallow_shp_dir, wallow_zips[i])
  download.file(paste0(geomac_wallow_prefix, wallow_zips[i]), 
                destfile)
  unzip(destfile, exdir = wallow_shp_dir)
  unlink(destfile)
  if (i == 1) {
    wallow_perims <- st_read(gsub('\\.zip', '.shp', destfile))
  } else {
    wallow_perims <- rbind(wallow_perims, 
                           st_read(gsub('\\.zip', '.shp', destfile)))
  }
}
wallow_perim_sf <- wallow_perims %>%
  arrange(-ACRES) %>%
  mutate(date = as.numeric(as.Date(DATE_, origin = '1960-10-01')), 
         date = date - min(date)) 

wallow_perim_plot <- wallow_perim_sf %>%
  ggplot(aes(fill = date)) +
  geom_sf(size = .1) + 
  scale_fill_viridis_c(option = 'A', breaks = c(1, 5, 10, 15, 20, 25), 
                       'Days since ignition') + 
  geom_sf(data = filter(wallow_perim_sf, DATE_ == as.Date('2011-05-31')), 
          fill = NA, color = 'white') + 
  theme_minimal() + 
  theme(panel.grid.minor = element_blank(), legend.position = 'none',
        axis.text.x = element_text(angle = 45)) + 
  annotate('text', x = -109.35, y = 33.7, label = 'May 31', color = 'white') + 
  xlab('') + 
  ylab('') + 
  ggtitle('A')
wallow_perim_plot

hectares_per_acre <- 0.404686
wallow_size_df <- wallow_perim_sf %>%
  as.data.frame %>%
  select(-geometry) %>%
  mutate(hectares = ACRES * hectares_per_acre, 
         date = DATE_, 
         value = hectares, 
         name = 'Burned area (hectares)', 
         pixel = 1, 
         label = 'Daily local conditions') %>%
  select(pixel, date, value, name, label)


pts_df <- daily_df %>%
  mutate(label = 'Daily local conditions') %>%
  full_join(wallow_size_df)

p2 <- pts_df %>%
  ggplot(aes(as.Date(date), value)) + 
  geom_point(aes(color = label), size = 0.6,
             data = filter(pts_df, name == 'Burned area (hectares)')) + 
  geom_point(aes(color = label), 
             size = .6, alpha = .02,
             data = filter(pts_df, name != 'Burned area (hectares)')) +
  geom_line(data = region_df, 
            aes(group = month, color = label), 
            size = 1) + 
  facet_wrap(~name, scales = 'free_y', ncol = 1) + 
  theme_minimal() + 
  xlab('') + 
  ylab('') + 
  theme(panel.grid.minor = element_blank(), legend.position = 'none') + 
  scale_x_date(date_labels = "%Y-%m-%d") + 
  scale_color_manual(values = c('black', region_color), 
                     guide = guide_legend(override.aes = list(
                       linetype = c("blank", "solid"),
                       shape = c(19, NA), 
                       alpha = c(1, 1))), 
                     "") + 
  ggtitle('B')
p2

wallow_perim_plot + p2 + plot_layout(widths = c(1, 1))
ggsave('fig/figure_12.pdf', width = 7, height = 5)
