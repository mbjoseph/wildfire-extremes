library(broom)
library(plotly)
library(ggthemes)

# Evaluating movement in 2D climate space for each ecoregion --------------
source('R/02-explore.R')

climate_d <- st_covs %>%
  select(NA_L3NAME, year, ym, cpr, crmin, ctmx, cvs, cpr12) %>%
  na.omit

pc <- climate_d %>%
  select(-NA_L3NAME, -year, -ym) %>%
  prcomp(scale = TRUE)

pc_d <- climate_d %>%
  bind_cols(as_tibble(pc$x)) %>%
  arrange(NA_L3NAME, ym) %>%
  left_join(er_df)

p <- pc_d %>%
  plot_ly(x = ~PC1, y = ~PC2, z = ~PC3, color = ~NA_L3NAME, text=~factor(ym),
             type = 'scatter3d', mode = 'lines',
             opacity = 1,
             line = list(width = 3, colorscale = 'Viridis')) %>%
  layout(showlegend = FALSE)
p

