# Loading and exploring some data
library(dplyr)
library(ggplot2)

d <- read.csv('data/d.csv', stringsAsFactors = FALSE)
# t: time
# y: observations


# Evaluate the distribution of y
d %>%
  ggplot(aes(x = y)) +
  geom_histogram(bins = 40)

d %>%
  ggplot(aes(x = y)) +
  facet_wrap(~t) +
  geom_histogram(bins = 40)

d %>%
  ggplot(aes(x = t, y = y)) +
  geom_jitter(alpha = .4) +
  theme_light()



# compute block maxima, grouping by t
extremes <- d %>%
  group_by(t) %>%
  summarize(max_y = max(y),
            min_y = min(y),
            n = n())

# visualize block maxima through time
extremes %>%
  ggplot(aes(x = t, y = max_y)) +
  geom_point() +
  xlab('t') +
  ylab('max(y)')

# is this due to more records in recent ts?
extremes %>%
  ggplot(aes(x = n, y = max_y)) +
  geom_point() +
  xlab('Number of records in data') +
  ylab('max(y)')

# if this is entirely driven by sampling,
# we might expect that smaller values are recorded in recent ts such that
# the min(y) decreases through time
extremes %>%
  ggplot(aes(x = t, y = min_y)) +
  geom_point() +
  xlab('t') +
  ylab('min(y)')
# nah, seems pretty stable

extremes %>%
  ggplot(aes(x = n, y = min_y)) +
  geom_point() +
  xlab('Number of records in data') +
  ylab('max(y)')
