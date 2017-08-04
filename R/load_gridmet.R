library(raster)
library(rgdal)
library(tidyverse)

gridmet_files <- list.files(file.path("data", "raw", "gridmet"),
                            pattern = "monthly",
                            recursive = TRUE,
                            full.names = TRUE)

# load a stack of all monthly gridmet layers
r <- stack(gridmet_files)

# convert to data frame for construction of model matrix
rd <- as.data.frame(r, xy = TRUE) %>%
  tbl_df %>%
  na.omit
