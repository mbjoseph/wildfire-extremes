library(sf)
library(tidyverse)
library(RColorBrewer)
library(patchwork)

mtbs <- read_rds('data/processed/mtbs.rds')
ecoregions <- read_rds('data/processed/ecoregions.rds')

l3_palette <- colorRampPalette(brewer.pal(12, 'Set3'))
l1_map <- ecoregions %>%
  group_by(NA_L1NAME) %>%
  summarize(n = n()) %>%
  ggplot() +
  geom_sf(aes(fill = NA_L1NAME), size = .1) + 
  theme_minimal() + 
  theme(panel.grid.major = element_line(colour = "lightgrey"), 
        legend.position = 'none') + 
  scale_fill_manual(values = l3_palette(length(unique(ecoregions$NA_L1NAME)))) + 
  ggtitle('B')

l2_map <- ecoregions %>% 
  group_by(NA_L2NAME) %>%
  summarize(n = n()) %>%
  ggplot() +
  geom_sf(aes(fill = NA_L2NAME), size = .1) + 
  theme_minimal() + 
  theme(panel.grid.major = element_line(colour = "lightgrey"), 
        legend.position = 'none') + 
  scale_fill_manual(values = l3_palette(length(unique(ecoregions$NA_L2NAME)))) + 
  ggtitle('C')

l3_map <- ecoregions %>%
  ggplot() +
  geom_sf(aes(fill = NA_L3NAME), size = .1) + 
  theme_minimal() + 
  theme(panel.grid.major = element_line(colour = "lightgrey"), 
        legend.position = 'none') + 
  scale_fill_manual(values = l3_palette(length(unique(ecoregions$NA_L3NAME)))) + 
  ggtitle('D')

fire_map <- ecoregions %>%
  ggplot() +
  geom_sf(size = .1, fill = 'white') + 
  geom_sf(data = mtbs, size = .2, inherit.aes = FALSE) + 
  theme_minimal() + 
  theme(panel.grid.major = element_line(colour = "lightgrey"), 
        legend.position = 'none') + 
  ggtitle('A')

p <- (fire_map  + l1_map) / (l2_map + l3_map)
ggsave(filename = 'fig/figure_2.pdf', plot = p, width = 10, height = 7)
