library(raster)
library(dplyr)

read_pdsi <- function(f) {
  pdsi <- stack(f, varname = "pdsisc")
  names(pdsi) <- paste(substr(names(pdsi), start = 2, stop = 5), 1:12, sep = "-")

  years_of_interest <- 1992:2013
  layers_to_extract <- expand.grid(y = paste0("X", years_of_interest),
                                   m = 1:12) %>%
    arrange(y) %>%
    apply(1, paste, collapse = ".", sep = "") %>%
    gsub(pattern = " ", replacement = "")
  pdsi <- subset(pdsi, layers_to_extract) %>%
    projectRaster(crs = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

  annual_means <- pdsi %>%
    stackApply(indices = rep(1:length(years_of_interest), each = 12), fun = mean)
  names(annual_means) <- years_of_interest
  annual_means
}

pdsi <- "data/raw/pdsisc.monthly.maps.1900-2099.r2.5x2.5.EnsAvg14Models.TP2.ipe=2.nc" %>%
  read_pdsi()
