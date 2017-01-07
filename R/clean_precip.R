library(tidyverse)
library(lubridate)
library(raster)

# precipitation data are available from NASA NEX-GDDP via planet os here:
# http://opennex.planetos.com/gddp/e8PE3

# to break up the massive us-precip.csv into multiple, use
# system("bash bash/parse-precip.sh")

# find processed precip files to parse
precip_files <- list.files(path = file.path("data", "processed"),
                           pattern = "precip-",
                           full.names = TRUE)
stopifnot(length(precip_files) > 0)

# extract column names from the big file
precip_headers <- readLines("data/raw/us-precip.csv", n = 1) %>%
  strsplit(split = ",") %>%
  unlist()

# loop over each year and compute annual means at each location
l <- list()
for (i in seq_along(precip_files)) {
  l[[i]] <- read_csv(file = precip_files[i], col_names = precip_headers) %>%
    group_by(Longitude, Latitude) %>%
    summarize(mean_precip = mean(Value),
              total_precip = sum(Value))
}
names(l) <- 1991:2013

# construct one big data frame and save it
precip_d <- bind_rows(l, .id = "year")
names(precip_d) <- tolower(names(precip_d))

precip_d %>%
  ungroup() %>%
  mutate(year = as.numeric(year)) %>%
  filter(year < 2014) %>%
  write_csv(path = "data/processed/cleaned-precip.csv")
