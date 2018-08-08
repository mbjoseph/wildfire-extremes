
# Plot a climate layer ----------------------------------------------------

library(raster)
library(viridis)
library(tidyverse)
source('R/helpers.R')

ecoregion_shp <- load_ecoregions()

r1 <- raster('data/processed/climate-data/rmin/monthly_rmin_2016.tif') %>%
  projectRaster(crs = crs(ecoregion_shp))

r2 <- raster('data/processed/climate-data/pr/monthly_pr_2016.tif') %>%
  projectRaster(crs = crs(ecoregion_shp))

r3 <- raster('data/processed/climate-data/vs/monthly_vs_2016.tif') %>%
  projectRaster(crs = crs(ecoregion_shp))

r4 <- raster('data/processed/climate-data/tmmx/monthly_tmmx_2016.tif') %>%
  projectRaster(crs = crs(ecoregion_shp))

d <- raster('data/processed/den10.tif') %>%
  projectRaster(to = r1) %>%
  mask(mask = r1)

values(d)[values(d) < 0] = NA

par(mfrow = c(5, 1), mar = rep(0, 4))
plot(r1, col = viridis(50), axes=FALSE, box=FALSE, legend = FALSE)
plot(log(r2), col = viridis(50), axes=FALSE, box=FALSE, legend = FALSE)
plot(r3, col = viridis(50), axes=FALSE, box=FALSE, legend = FALSE)
plot(r4, col = viridis(50), axes=FALSE, box=FALSE, legend = FALSE)
plot(log(d), col = viridis(50, option = 'C'), axes = FALSE, box = FALSE, legend = FALSE)

