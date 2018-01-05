
# Helper functions for fire extremes project ------------------------------

load_ecoregions <- function() {
  if (!dir.exists("data/raw/us_eco_l3/")) {
    download.file("ftp://newftp.epa.gov/EPADataCommons/ORD/Ecoregions/us/us_eco_l3.zip",
                  destfile = "data/raw/us_eco_l3.zip")
    unzip("data/raw/us_eco_l3.zip",
          exdir = "data/raw/us_eco_l3/")
  }

  ecoregion_shp <- rgdal::readOGR(dsn = "data/raw/us_eco_l3/",
                                   layer = "us_eco_l3")

  ecoregion_shp
}
