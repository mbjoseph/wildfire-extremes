library(raster)
library(rgdal)
library(tidyverse)

maca_files <- list.files("data/raw/maca", pattern = ".nc", full.names = TRUE)
maca_types <- c("huss", "pr", "rsds", "tasmax", "tasmin", "was")

rlist <- list()
k <- 1 # list index

for (i in seq_along(maca_types)) {
  which_files <- maca_files[grepl(pattern = maca_types[i], x = maca_files)]

  for (j in seq_along(which_files)) {
    rlist[[k]] <- stack(which_files[j])
    names(rlist[[k]]) <- gsub(pattern = "X",
                              replacement = paste0(maca_types[i], "_"),
                              x = names(rlist[[k]]))
    k <- k + 1
  }
}

r <- stack(rlist)


# subset to only keep relevant years
to_exclude <- c(1970:1983, 2016:2025)
for (i in seq_along(to_exclude)) {
  keep <- grep(pattern = to_exclude[i], x = names(r), invert = TRUE, value = TRUE)
  r <- subset(r, keep)
}

# generate annual mean values
r_year <- strsplit(names(r), split = "_") %>%
  lapply(`[[`, 2) %>%
  unlist() %>%
  substr(1, 4) %>%
  parse_number

variable <- strsplit(names(r), split = "_") %>%
  lapply(`[[`, 1) %>%
  unlist()

annual_r <- list()
k <- 1
for (i in seq_along(unique(r_year))) {
  for (j in seq_along(unique(variable))) {
    which_idx <- r_year == unique(r_year)[i] &
      variable == unique(variable)[j]
    annual_r[[k]] <- subset(r, names(r)[which_idx]) %>%
      mean
    k <- k + 1
  }
}

annual_rstack <- stack(annual_r)
names(annual_rstack) <- expand.grid(unique(variable), unique(r_year)) %>%
  as.matrix() %>%
  apply(1, paste, collapse = "_")

coarse_r <- aggregate(annual_rstack[[1]], fact = 4)
poly_grid <- rasterToPolygons(coarse_r)

