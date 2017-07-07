library(raster)
library(lubridate)
library(rgdal)
library(tidyverse)

usa_shp <- readOGR(dsn = "data/raw/cb_2016_us_nation_20m",
                   layer = "cb_2016_us_nation_20m")

usa_shp <- spTransform(usa_shp,
                       CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))

summarize_by_month <- function(file, mask_shp) {
  # generate monthly summaries from daily climate data
  file_split <- file %>%
    basename %>%
    strsplit(split = "_") %>%
    unlist
  var <- file_split[1]
  year <- substr(file_split[2], start = 1, stop = 4)

  out_name <- gsub(file,
                   pattern = basename(file),
                   replacement = paste0("monthly_",
                                        basename(file))) %>%
    gsub(pattern = ".nc", replacement = ".tif")
  if (file.exists(out_name)) {
    return(paste("File", out_name, "already exists"))
  }
  if (as.numeric(year) < 1983) {
    return("Year outside range of consideration")
  }

  # determine which function to use
  fun <- ifelse(var == "pr", sum, mean)

  # otherwise, generate monthly summary
  raster <- stack(file)
  start_date <- as.Date(paste(year, "01", "01", sep = "-"))
  end_date <- as.Date(paste(year, "12", "31", sep = "-"))
  date_seq <- seq(start_date, end_date, by = "1 day")
  month_seq <- month(date_seq)

  res <- stackApply(raster, month_seq, fun = fun)
  corrected_res <- flip(t(res), direction = "x")
  names(corrected_res) <- paste(var, year,
                                unique(month(date_seq, label = TRUE)),
                                sep = "_")
  p4string <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  projection(corrected_res) <- CRS(p4string)
  masked_res <- mask(corrected_res, mask_shp)
  writeRaster(masked_res, out_name, format = "GTiff")
  return(paste("File", out_name, "written"))
}

# identify which files need to be processed
daily_files <- list.files("data/raw/gridmet", pattern = "nc",
                          recursive = TRUE, full.names = TRUE)

system.time(res <- summarize_by_month(daily_files[[7]], usa_shp))

res

system.time(monthly_precip <- monthly_summary("pr", "1983", fun = sum, usa_shp))

plot(monthly_precip)

