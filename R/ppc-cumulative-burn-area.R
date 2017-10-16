library(tidyverse)
# Posterior predictive check for cumulative burn area ---------------------

mtbs_hist <- read_rds('mtbs_hist.rds')
hist_df <- read_rds('ppc_hist.rds')


mtbs_hist_df <- tibble(count = mtbs_hist$counts,
                       midpoint = mtbs_hist$mids) %>%
  mutate(approx_burn_area = count * midpoint,
         cum_burn_area = cumsum(approx_burn_area),
         iter = 1)


hist_df %>%
  ggplot(aes(log(midpoint), log(cum_burn_area), group = iter)) +
  geom_line(alpha = .1) +
  geom_line(data = mtbs_hist_df, color = 'red') +
  coord_cartesian(ylim = c(16.5, 25),
                  xlim = c(7.5, 20)) +
  xlab('log(Burn area)') +
  ylab('log(Cumulative burn area)')

