
# Fetch historical CMIP6 land cover data ----------------------------------
# from http://luh.umd.edu/data.shtml
url <- "https://www.dropbox.com/sh/eanwzwapcfpue1s/AAB2hjnmQ4UFlmmn6b_Btioxa/states.nc?dl=1"
dest <- "data/raw/cmip6/states.nc"
if (!file.exists(dest)) {
  download.file(url, destfile = dest)
}
