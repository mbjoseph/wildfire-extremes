library(assertthat)

# Acquiring raw data ------------------------------------------------------
# This script will check to see whether all of the raw data files exist locally
# and download any missing files.

raw_prefix <- file.path("data", "raw")

# Level 3 ecoregion data from EPA -----------------------------------------

ecoregion_prefix <- file.path(raw_prefix, "us_eco_l3")
ecoregion_shp <- file.path(ecoregion_prefix, "us_eco_l3.shp")

if (!file.exists(ecoregion_shp)) {
  loc <- "ftp://newftp.epa.gov/EPADataCommons/ORD/Ecoregions/us/us_eco_l3.zip"
  dest <- paste0(ecoregion_prefix, ".zip")
  download.file(loc, dest)
  unzip(dest, exdir = ecoregion_prefix)
  unlink(dest)
  assert_that(file.exists(ecoregion_shp))
}



# MTBS fire data ----------------------------------------------------------

mtbs_prefix <- file.path(raw_prefix, "mtbs_fod_pts_data")
mtbs_shp <- file.path(mtbs_prefix, 'mtbs_fod_pts_20170501.shp')

if (!file.exists(mtbs_shp)) {
  loc <- "https://edcintl.cr.usgs.gov/downloads/sciweb1/shared/MTBS_Fire/data/composite_data/fod_pt_shapefile/mtbs_fod_pts_data.zip"
  dest <- paste0(mtbs_prefix, ".zip")
  download.file(loc, dest)
  unzip(dest, exdir = mtbs_prefix)
  unlink(dest)
  assert_that(file.exists(mtbs_shp))
}
