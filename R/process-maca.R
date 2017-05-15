library(raster)
library(rgdal)

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

writeRaster(r, "data/processed/maca.tif")
