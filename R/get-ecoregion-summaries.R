library(raster)
library(snowfall)
library(rgdal)

if (!dir.exists("data/raw/us_eco_l3/")) {
  download.file("ftp://newftp.epa.gov/EPADataCommons/ORD/Ecoregions/us/us_eco_l3.zip",
                destfile = "data/raw/us_eco_l3.zip")
  unzip("data/raw/us_eco_l3.zip",
        exdir = "data/raw/us_eco_l3/")
}

ecoregion_shp <- readOGR(dsn = "data/raw/us_eco_l3/",
                   layer = "us_eco_l3")

ecoregion_shp <- spTransform(ecoregion_shp,
                       CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))


tifs <- list.files(pattern = ".tif")

extract_one <- function(filename, ecoregion_shp) {
  raster::extract(raster::stack(filename), ecoregion_shp,
                  na.rm = TRUE, fun = mean, df = TRUE)
}

sfInit(parallel = TRUE, cpus = 4)
sfExport(list = c("ecoregion_shp"))

# try on subset
extractions <- sfLapply(as.list(tifs),
                        fun = extract_one,
                        ecoregion_shp = ecoregion_shp)
sfStop()
