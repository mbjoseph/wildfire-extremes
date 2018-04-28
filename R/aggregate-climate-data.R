library(raster)
library(lubridate)
library(rgdal)
library(tidyverse)
library(parallel)
library(pbapply)

if (!dir.exists("data/raw/cb_2016_us_nation_20m/")) {
  download.file("http://www2.census.gov/geo/tiger/GENZ2016/shp/cb_2016_us_nation_20m.zip",
                destfile = "data/raw/cb_2016_us_nation_20m.zip")
  unzip("data/raw/cb_2016_us_nation_20m.zip",
        exdir = "data/raw/cb_2016_us_nation_20m/")
}

usa_shp <- readOGR(dsn = "data/raw/cb_2016_us_nation_20m",
                   layer = "cb_2016_us_nation_20m")

usa_shp <- spTransform(usa_shp,
                       CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))

summarize_by_month <- function(url, mask_shp) {
  # takes a climate data url and a masking spatial polygon as input
  # and makes a monthly summary of the climate data, masked by the polygon


  # pull the raw climate data file down locally
  file <- basename(url)

  # generate an output filename
  nc_name <- gsub(file,
                  pattern = basename(file),
                  replacement = paste0("monthly_",
                                       basename(file)))
  out_name <- gsub(x = nc_name, pattern = ".nc", replacement = ".tif")

  # determine which climate variable, and which year we have
  file_split <- unlist(strsplit(basename(file), split = "_"))
  var <- file_split[1]
  year <- substr(file_split[2], start = 1, stop = 4)

  # pass over file if necessary
  if (basename(out_name) %in% list.files(pattern = basename(out_name),
                                         recursive = TRUE)) {
    return(paste("File", out_name, "already exists"))
  }
  if (as.numeric(year) < 1983) {
    return("Year outside range of consideration")
  }

  # if we haven't already exited, download the original data file
  download.file(url = url, destfile = file)

  # determine which function to use to aggregate
  if (var == "pr") {
    fun <- sum
  } else {
    fun <- mean
  }

  # generate monthly summary
  raster <- raster::stack(file)
  start_date <- as.Date(paste(year, "01", "01", sep = "-"))
  end_date <- as.Date(paste(year, "12", "31", sep = "-"))
  date_seq <- seq(start_date, end_date, by = "1 day")
  date_seq <- date_seq[1:raster::nlayers(raster)]
  month_seq <- lubridate::month(date_seq)

  res <- raster::stackApply(raster, month_seq, fun = fun)
  corrected_res <- raster::flip(raster::t(res), direction = "x")
  names(corrected_res) <- paste(var, year,
                                unique(lubridate::month(date_seq, label = TRUE)),
                                sep = "_")
  p4string <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  raster::projection(corrected_res) <- sp::CRS(p4string)
  masked_res <- raster::mask(corrected_res, mask_shp)
  raster::writeRaster(masked_res, out_name, format = "GTiff")
  unlink(file)
  return(paste("File", out_name, "written"))
}

climate_data_urls <- read.csv('data/raw/climate-data.csv',
                              stringsAsFactors = FALSE)


# Summarize daily climate data by month ---------------------------
print('Aggregating daily climate data to monthly means. May take a while...')
pboptions(type = 'txt', use_lb = TRUE)
cl <- makeCluster(getOption("cl.cores", detectCores() / 2))
pblapply(X = climate_data_urls$url, 
         FUN = summarize_by_month, 
         mask_shp = usa_shp, 
         cl = cl)
stopCluster(cl)

# Move files to proper directories -------------------------------------------
tifs_in_home_dir <- list.files(pattern = ".tif") %>%
  tibble(filename = .) %>%
  separate(filename, into = c("time_interval", "variable", "year"), sep = "_") %>%
  mutate(current_name = paste(time_interval, variable, year, sep = "_"))

variables <- distinct(tifs_in_home_dir, variable) %>%
  unlist()

# create dirs for each variable
dirs_to_make <- file.path("data", "processed", "climate-data", variables)
sapply(dirs_to_make, dir.create, recursive = TRUE)

# move each tif file to the proper location
tifs_in_home_dir %>%
  mutate(dest_dir = file.path("data", "processed", 'climate-data', variable),
         dest_file = file.path(dest_dir, current_name)) %>%
  do(file.rename(.$current_name, .$dest_file) %>% tibble)

print('Daily climate data aggregated to monthly summaries.')
