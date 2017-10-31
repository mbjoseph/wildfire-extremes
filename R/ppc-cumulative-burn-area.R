library(tidyverse)
# Posterior predictive check for cumulative burn area ---------------------

mtbs_hist <- read_rds('mtbs_hist.rds')
hist_df <- read_rds('ppc_hist.rds')


mtbs_hist_df <- tibble(count = mtbs_hist$counts,
                       midpoint = mtbs_hist$mids) %>%
  mutate(approx_burn_area = count * midpoint,
         cum_burn_area = cumsum(approx_burn_area),
         iter = 1,
         p_burn_area = cum_burn_area / max(cum_burn_area))


hist_df %>%
  ggplot(aes(log(midpoint), log(cum_burn_area), group = iter)) +
  geom_line(alpha = .1) +
  geom_line(data = mtbs_hist_df, color = 'red') +
  coord_cartesian(ylim = c(10.5, 21),
                  xlim = c(5, 16)) +
  xlab('log(Burn area)') +
  ylab('log(Cumulative burn area)')

# normalized by total burn area
hist_df %>%
  group_by(iter) %>%
  mutate(total_burn_area = max(cum_burn_area),
         p_burn_area = cum_burn_area / total_burn_area) %>%
  ggplot(aes(log(midpoint), p_burn_area, group = iter)) +
  geom_line(alpha = .1) +
  geom_line(data = mtbs_hist_df, color = 'red') +
  # coord_cartesian(ylim = c(10.5, 21),
  #                 xlim = c(5, 16)) +
  xlab('log(Burn area)') +
  scale_x_log10() +
  ylab('Cumulative proportion of total burn area')
