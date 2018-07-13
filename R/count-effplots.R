
# Count plots -------------------------------------------------------------

library(tidyverse)
library(rstan)
library(ggthemes)
library(ggrepel)
library(sf)
library(Matrix)

colnamesX <- read_rds('data/processed/colnamesX.rds')
st_covs <- read_rds('data/processed/st_covs.rds')
vars <- read_rds('data/processed/vars.rds')
X <- read_rds('data/processed/X.rds')
ecoregion_df <- read_rds('data/processed/ecoregions.rds') %>%
  as("Spatial") %>%
  data.frame

# Extract posterior draws and visualize results ---------------------------

post <- rstan::extract(read_rds(path = 'zinb_full_fit.rds'), 
                       pars = c('Rho_beta', 'beta', 'lp__'))


# Coefficient caterpillar plot ------------------------------

beta_summary <- post$beta %>%
  reshape2::melt(varnames = c('iter', 'dim', 'col')) %>%
  tbl_df %>%
  group_by(dim, col) %>%
  summarize(median = median(value),
            mean = mean(value),
            lo = quantile(value, 0.025),
            hi = quantile(value, 0.975),
            p_neg = mean(value < 0),
            p_pos = mean(value > 0)) %>%
  ungroup %>%
  mutate(variable = colnamesX[col],
         nonzero = p_neg > .87 | p_pos > .87, 
         Response = factor(dim,
                           labels = c(expression(paste('Negative binomial coefficients: ', beta^(mu))),
                                      expression(paste('Zero-inflation coefficients: ', beta^(pi))))),
         variable = gsub('bs_', '', variable),
         variable = gsub('rmin', 'humidity', variable),
         variable = gsub('tmmx', 'temperature', variable),
         variable = gsub('vs', 'wind speed', variable),
         variable = gsub('prev_12mo_precip', '12 mo. precip.', variable),
         variable = gsub('^pr', 'precipitation', variable),
         variable = gsub('log_housing_density', 'Housing density', variable),
         variable = gsub(':', ' x ', variable),
         variable = gsub('NA_', '', variable),
         variable = gsub('NAME', ' ', variable),
         variable = gsub('_', ': ', variable),
         variable = gsub('USA', '', variable),
         variable = gsub('AND', '&', variable),
         variable = tolower(variable),
         variable = tools::toTitleCase(variable))

# show all coefficients
coef_d <- beta_summary %>%
  group_by(Response) %>%
  mutate(max_abs = max(abs(median)),
         rel_size = median / max_abs)

coefplot <- coef_d %>%
  ggplot(aes(x = median,
             y = variable,
             color = ifelse(nonzero, sign(median), 0))) +
  geom_point(aes(size = nonzero), alpha = .6) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major.y = element_blank()) +
  facet_wrap(~ Response, scales = 'free', nrow = 1, labeller = label_parsed) +
  geom_segment(aes(x = lo, xend = hi,
                   yend = variable), 
               size = .5, alpha = .1,
               data = filter(coef_d, !nonzero)) +
  geom_segment(aes(x = lo, xend = hi,
                   yend = variable), 
               alpha = .6,
               data = filter(coef_d, nonzero)) +
  scale_size_manual(values = c(.2, 1)) +
  ylab('') +
  xlab('Coefficient value') +
  scale_color_gradient2(mid = 'lightgrey', low = 'dodgerblue', high = 'tomato') +
  theme(legend.position = 'none') +
  geom_text_repel(aes(label = ifelse(nonzero, variable, '')),
                  color = 'black', size = 2)
ggsave('fig/all-coefs.png', plot = coefplot, width = 6, height = 3.5)


# Partial effect plots ----------------------------------------------------

partial_effs <- list()
n_iter <- length(post$lp__)
unique_ers <- unique(st_covs$NA_L3NAME)
pb <- txtProgressBar(max = length(unique_ers), style = 3)
for (v in seq_along(vars)) {
  print(paste('Processing', vars[v]))
  for (i in seq_along(unique_ers)) {
    setTxtProgressBar(pb, i)
    X_sub <- X[st_covs$NA_L3NAME == unique_ers[i], ]
    cols <- grepl(vars[v], colnames(X_sub))
    if (vars[v] == 'pr') {
      cols <- grepl(vars[v], colnames(X_sub)) &
        !grepl('12mo', colnames(X_sub))
    }

    effects <- X_sub[, cols] %*% t(post$beta[, 1, cols]) %>% # 1 --> negbinom mean
      apply(1, quantile, c(0.025, .5, .975))
    
    subdf <- st_covs %>%
      filter(NA_L3NAME == unique_ers[i]) %>%
      mutate(row_id = parse_number(rownames(X_sub)))
    
    partial_effs[[length(partial_effs) + 1]] <- effects %>%
      reshape2::melt(varnames = c('quantity', 'row_id')) %>%
      as_tibble %>%
      mutate(quantity = case_when(.$quantity == "2.5%" ~ 'lo',
                                  .$quantity == "50%" ~ 'med',
                                  .$quantity == "97.5%" ~ 'hi'),
             var = vars[v]) %>%
      spread(quantity, value) %>%
      left_join(select(subdf, row_id, NA_L3NAME, ym))
  }
}
close(pb)


name_df <- tibble(covariate = vars,
                  fancy_name = c('Housing density log(units/km^2)',
                                 'Wind speed (m/s)',
                                 'Precipitation: same month (mm)',
                                 'Precipitation: 12 month (mm)',
                                 'Air temperature (C)',
                                 'Relative humidity (%)'))

effect_plot_df <- list()
for (i in seq_along(vars)) {
  effect_plot_df[[i]] <- partial_effs %>%
    bind_rows %>%
    filter(var == vars[i]) %>%
    left_join(st_covs) %>%
    mutate(covariate = vars[i]) %>%
    left_join(name_df) %>%
    rename(covariate_value = !!vars[i]) %>%
    select(covariate, covariate_value, med, lo, hi,
           NA_L3NAME, ym, NA_L1NAME, NA_L2NAME,
           fancy_name)
}

better_name <- distinct(ecoregion_df, NA_L1NAME) %>%
  mutate(l1_er = tools::toTitleCase(tolower(NA_L1NAME))) %>%
  mutate(ifelse(l1_er == 'Southern Semi-Arid Highlands',
                'Southern Semi-arid Highlands',
                l1_er))


regions_to_highlight <- c("Eastern Temperate Forests", 
                          'Great Plains', 
                          'North American Deserts')

plot_data <- effect_plot_df %>%
  bind_rows %>% 
  mutate(covariate_value = case_when(
    .$covariate == 'tmmx' ~ covariate_value - 273.15,
    TRUE ~ covariate_value)) %>%
  left_join(better_name) %>%
  mutate(l1_er_subset = ifelse(l1_er %in% regions_to_highlight, 
                               l1_er, "Other"))



pal <- c("#D55E00", "#56B4E9", "#009E73", scales::alpha(1, .2))#"#999999")

p <- plot_data %>%
  ggplot(aes(covariate_value, y = med)) +
  theme_minimal() +
  geom_line(aes(group = NA_L3NAME, color = l1_er_subset)) +
  ylab('Partial effect') +
  facet_wrap(~fancy_name, scales = 'free_x', ncol = 3, 
             strip.position = 'bottom') +
  scale_color_manual(values = pal, '') +
  xlab('') +
  theme(panel.grid.minor = element_blank(),
  legend.position = 'bottom') + 
  scale_y_continuous(breaks = c(-3, 0, 3))
p
ggsave('fig/count-partial-effs.png', plot = p, width = 6, height = 4)



# Save numeric summaries for text -----------------------------------------

post$Rho_beta[, 1, 2] %>%
  as_tibble %>%
  write_csv('data/processed/rho_beta.csv')


