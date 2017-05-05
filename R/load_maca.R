library(raster)
library(dplyr)

# Loading MACA climate data -----------------------------------------------
ncdf_files <- list.files("data/raw/maca", pattern = ".nc", full.names = TRUE)

wind_speed <- ncdf_files %>%
  grep(pattern = "was", value = TRUE) %>%
  grep(pattern = "rcp45", value = TRUE, invert = TRUE) %>%
  stack
