library(cmipr)
library(tidyverse)
library(raster)

files <- cmip_list_files('bcsd/yearly/ncar_ccsm3_0.7')

files <- files %>%
  filter(grepl("sresa1b.monthly.Prcp", file))

res <- lapply(file.path('bcsd/yearly/ncar_ccsm3_0.7', files$file),
              cmip_fetch)

prcp_layers <- lapply(res, FUN = grepl,
                      pattern = "sresa1b\\.monthly\\.Prcp") %>%
  unlist %>%
  which

rlist <- res %>%
  lapply(FUN = cmip_read)

rstack <- raster::stack(rlist)

names(rstack) <- files$file
