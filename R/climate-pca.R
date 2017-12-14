library(broom)
library(plotly)
library(ggthemes)

# Evaluating movement in 2D climate space for each ecoregion --------------
source('R/02-explore.R')

climate_d <- st_covs %>%
  select(NA_L3NAME, year, ym, cpr, crmin, ctmx, cvs, cpr12) %>%
  gather(variable, value, -NA_L3NAME, -year, -ym) %>%
  mutate(month = substr(ym, 1, 3)) %>%
  unite(variable, variable, month) %>%
  select(-ym) %>%
  spread(variable, value) %>%
  na.omit

pc <- climate_d %>%
  select(-NA_L3NAME, -year) %>%
  prcomp(scale = TRUE)

pc_d <- climate_d %>%
  bind_cols(as_tibble(pc$x)) %>%
  arrange(NA_L3NAME, year) %>%
  left_join(er_df)

p <- pc_d %>%
  ggplot(aes(PC1, PC2, color = NA_L3NAME, group = NA_L3NAME)) +
  geom_path(alpha = .8) +
  geom_point(aes(group = year), alpha = .5, size = .1) +
  theme_minimal() +
  theme(legend.position = 'none')

ggplotly(p)


p <- plot_ly(pc_d, x = ~PC1, y = ~PC2, z = ~PC3, color = ~NA_L3NAME,
             type = 'scatter3d', mode = 'lines',
             opacity = 1,
             line = list(width = 3, colorscale = 'Viridis')) %>%
  layout(showlegend = FALSE)
p


p <- plot_ly(pc_d, x = ~PC1, y = ~PC2, z = ~PC3, color = ~NA_L3NAME,
             sizes = 2) %>%
  add_lines()
p
