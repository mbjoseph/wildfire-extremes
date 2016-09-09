
# fit threshold models for each l3 ecoregion ------------------------------
source('R/thresholds/clean_data.R')
library(evmix)
hist(d$fire_size, breaks = 100)

d <- subset(d, stat_cau_1 == 'Lightning')

fits <- vector(mode = "list", length = length(unique(d$us_l3name)))
unique_regions <- unique(d$us_l3name)
thresh_df <- data.frame(us_l3name = unique_regions,
                        threshold = NA, xi = NA, n = NA, conv = NA)

for (i in seq_along(unique_regions)) {
  subd <- subset(d, us_l3name == unique_regions[i])
  if (nrow(subd) > 20) {
    fits[[i]] <- fweibullgpd(subd$fire_size)
    thresh_df$conv[i] <- fits[[i]]$conv
    if (fits[[i]]$conv) {
      thresh_df$threshold[i] <- fits[[i]]$u
      thresh_df$xi[i] <- fits[[i]]$xi
    }
    thresh_df$n[i] <- nrow(subd)
  }
}

thresh_df %>%
  arrange(-threshold)

## visualize thresholds -------------------------------
thresh_df %>%
  ggplot(aes(x = threshold)) +
  geom_rug(color = 'red') +
  geom_density(aes(x = fire_size, group = us_l3name), data = d) +
  scale_x_log10() +
  theme_minimal() +
  xlab('Fire size') +
  ylab('Fire size density')

thresh_df %>%
  ggplot(aes(x = n, y = threshold)) +
  geom_point() +
  xlab('Number of fires in ecoregion') +
  ylab('Estimated threshold')

thresh_df %>%
  ggplot(aes(x = n, y = xi)) +
  geom_point() +
  geom_abline(slope = 0, intercept = 0, linetype = 'dashed') +
  xlab('Number of fires in ecoregion') +
  ylab('Estimated GPD shape parameter')

thresh_df %>%
  ggplot(aes(x = threshold, y = xi)) +
  geom_point() +
  geom_abline(slope = 0, intercept = 0, linetype = 'dashed') +
  xlab('Estimated threshold') +
  ylab('Estimated GPD shape parameter')



## Compute and plot estimated pdfs ------------------------------------------
lo <- 1000
xseq <- exp(seq(min(log(d$fire_size)), log(100000), length.out = lo))
pdf_df <- data.frame(fire_size = rep(xseq, times = length(unique_regions)),
                     p_fire_size = NA,
                     us_l3name = rep(unique_regions, each = lo))

for (i in 1:nrow(thresh_df)) {
  if (!is.na(thresh_df$n[i])) {
    indx <- pdf_df$us_l3name == unique_regions[i]
    pdf_df$p_fire_size[indx] <- dweibullgpd(xseq,
                                         wshape = fits[[i]]$wshape,
                                         wscale = fits[[i]]$wscale,
                                         u = fits[[i]]$u,
                                         xi = fits[[i]]$xi)
  }
}

pdf_df %>%
  filter(!is.na(p_fire_size)) %>%
  ggplot() +
  geom_line(aes(x = fire_size, y = p_fire_size, group = us_l3name)) +
  scale_x_log10() +
#  scale_y_log10() +
  geom_jitter(aes(x = fire_size, y = -0.002), data = d, alpha = .3, size = .5,
              width = 0, height = .001) +
  geom_rug(aes(x = threshold), data = thresh_df, col = 'red')

# try on log scale to see upper bounds when xi < 0

## map it! -----------------------------------------------------
theme_map <- theme(axis.line = element_blank(),
                   axis.text.x = element_blank(),
                   axis.text.y = element_blank(),
                   axis.ticks = element_blank(),
                   axis.title.x = element_blank(),
                   axis.title.y = element_blank(),
                   panel.background = element_blank(),
                   panel.border = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   plot.background = element_blank())

er_df %>%
  full_join(thresh_df) %>%
  ggplot(aes(long, lat, group = group)) +
    geom_polygon(aes(fill = threshold)) +
    geom_polygon(color = 'black', alpha = .1, size = .1) +
    coord_equal() +
    scale_fill_gradientn(colors = c('midnightblue', 'darkblue',
                                    'blue', 'dodgerblue', 'cyan'),
                         'Estimated threshold') +
    theme_map

er_df %>%
  full_join(thresh_df) %>%
  ggplot(aes(long, lat, group = group)) +
  geom_polygon(aes(fill = xi)) +
  geom_polygon(color = 'black', alpha = .1, size = .1) +
  coord_equal() +
  scale_fill_gradientn(colors = c('black', 'darkred', 'red3', 'red',
                                  'orange', 'yellow'),
                       'Estimated shape parameter') +
  theme_map

# now plot the observed number of extreme fires in each ecoregion
full_join(thresh_df, d) %>%
  mutate(is_extreme = fire_size > threshold) %>%
  group_by(us_l3name) %>%
  summarize(n_extremes = sum(is_extreme),
            extreme_density = sum(is_extreme) / sum(shape_area)) %>%
  full_join(er_df) %>%
  ggplot(aes(long, lat, group = group)) +
  geom_polygon(aes(fill = n_extremes)) +
  #geom_polygon(aes(fill = extreme_density)) +
  geom_polygon(color = 'black', alpha = .1, size = .1) +
  coord_equal() +
  scale_fill_gradientn(colors = c('black', 'darkred', 'red3', 'red',
                                  'orange', 'yellow'),
                       'Extremely large fires\n1992-2013') +
  theme_map

overall_threshold <- mean(thresh_df$threshold, na.rm = TRUE)
