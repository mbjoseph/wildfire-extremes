
# Count plots -------------------------------------------------------------

source('R/02-explore.R')
library(extraDistr)
library(ggridges)
library(ggthemes)

fit <- read_rds(path = list.files(pattern = 'nb_fit.*'))

# Evaluate convergence ----------------------------------------------------
traceplot(fit, inc_warmup = TRUE)

traceplot(fit, pars = c('tau', 'c', 'alpha'))
traceplot(fit, pars = c('Rho_beta'))

# Extract posterior draws and visualize results ---------------------------
post <- rstan::extract(fit)
str(post)
rm(fit)
gc()



# Coefficients that seem to be far from zero ------------------------------
beta_df <- post$beta %>%
  reshape2::melt(varnames = c('iter', 'dim', 'col')) %>%
  tbl_df

beta_summary <- beta_df %>%
  group_by(dim, col) %>%
  summarize(median = median(value),
            lo = quantile(value, 0.025),
            hi = quantile(value, 0.975),
            p_neg = mean(value < 0),
            p_pos = mean(value > 0)) %>%
  ungroup %>%
  mutate(variable = colnamesX[col],
         nonzero = p_neg > .9 | p_pos > .9)

beta_summary %>%
  filter(nonzero) %>%
  select(dim, col, p_neg, p_pos, variable) %>%
  left_join(beta_df) %>%
  ggplot(aes(value, variable)) +
  geom_density_ridges() +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  facet_wrap(~ dim, scales = 'free') +
  theme(axis.text.y = element_text(size = 7)) +
  ggtitle('Estimated effects on # of fires')
ggsave('fig/fire-effs.png', width = 8, height = 8)

rm(beta_df)
gc()

