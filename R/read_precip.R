library(tidyverse)
library(viridis)
library(lubridate)

d <- read_csv("~/Downloads/data.csv") %>%
  mutate(year = year(Date)) %>%
  group_by(Longitude, Latitude, year) %>%
  summarize(value = mean(Value))

d %>%
  ggplot(aes(x = Longitude, y = Latitude, fill = value)) +
  geom_raster() +
  scale_fill_viridis() +
  facet_wrap(~ year)
