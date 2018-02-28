
# Count plots -------------------------------------------------------------

source('R/02-explore.R')
source('R/make-stan-d.R')
library(ggridges)
library(viridis)
library(ggthemes)
library(plotly)
library(ggrepel)

fit <- read_rds(path = 'zinb_full_fit.rds')

# Evaluate convergence ----------------------------------------------------
traceplot(fit, inc_warmup = TRUE)

traceplot(fit, pars = c('tau', 'c', 'alpha', 'Rho_beta'))

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
            mean = mean(value),
            lo = quantile(value, 0.025),
            hi = quantile(value, 0.975),
            p_neg = mean(value < 0),
            p_pos = mean(value > 0)) %>%
  ungroup %>%
  mutate(variable = colnamesX[col],
         nonzero = p_neg > .85 | p_pos > .85) %>%
  mutate(Response = ifelse(dim == 1,
                           'Negative binomial component',
                           'Zero-inflation component'),
         variable = gsub('bs_', '', variable),
         variable = gsub('crmin', 'humidity', variable),
         variable = gsub('ctmx', 'temperature', variable),
         variable = gsub('cvs', 'wind speed', variable),
         variable = gsub('cpr12', '12 mo. precip.', variable),
         variable = gsub('cpr', 'precipitation', variable),
         variable = gsub('chd', 'Housing density', variable),
         # variable = ifelse(grepl(':', x = variable),
         #                   paste0('(', variable, ')'),
         #                   variable),
         variable = gsub(':', ' x ', variable),
         variable = gsub('NA_', '', variable),
         variable = gsub('NAME', ' ', variable),
         variable = gsub('_', ': ', variable),
         variable = tolower(variable),
         variable = tools::toTitleCase(variable))

# show all coefficients
coefplot <- beta_summary %>%
  group_by(Response) %>%
  mutate(max_abs = max(abs(median)),
         rel_size = median / max_abs) %>%
  ggplot(aes(y = median,
             x = variable,
             color = ifelse(nonzero, sign(median), rel_size))) +
  geom_point(aes(size = nonzero)) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank()) +
  facet_wrap(~ Response, scales = 'free', nrow = 2) +
  geom_segment(aes(y = lo, yend = hi,
                   xend = variable, size = nonzero)) +
  scale_size_manual(values = c(.2, 1)) +
  xlab('Coefficient') +
  ylab('Coefficient value') +
  scale_color_gradient2(mid = 'lightgrey', low = 'blue', high = 'red') +
  theme(legend.position = 'none') +
  geom_text_repel(aes(label = ifelse(nonzero, variable, '')),
                  alpha = 1, size = 2.5, color = 'black')
coefplot
ggsave('fig/all-coefs.png', width = 8, height = 4)

# bivariate effects plot
bivar_d <- beta_summary %>%
  select(-dim, -col, -p_neg, -p_pos) %>%
  gather(which_val, value, -Response, -variable) %>%
  mutate(response_var = paste(Response, which_val, sep = '_')) %>%
  select(-Response, -which_val) %>%
  spread(response_var, value)

ptsize <- .5
bivar_effs <- bivar_d %>%
  ggplot(aes(`Zero-inflation component_median`, `Negative binomial component_median`)) +
  geom_point(alpha = .1, col = 'grey', shape = 1, size = ptsize) +
  geom_point(data = filter(bivar_d, `Negative binomial component_nonzero` |
                                      `Zero-inflation component_nonzero`),
             size = ptsize) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank()) +
  geom_text_repel(aes(label = variable),
                  data = filter(bivar_d, `Negative binomial component_nonzero` |
                                  `Zero-inflation component_nonzero`),
                  size = 3, max.iter = 10000, segment.alpha = .2) +
  ylab('Negative binomial effect') +
  xlab('Zero-inflation effect')
bivar_effs
ggsave('fig/bivar-effs.png', width = 6, height = 4)

# show important coefficients
beta_summary %>%
  filter(nonzero) %>%
  select(dim, col, p_neg, p_pos, variable,
         median, mean, variable, Response) %>%
  left_join(beta_df) %>%
  ggplot(aes(value, reorder(variable, mean), fill = median)) +
  theme_minimal() +
  geom_density_ridges(scale = 3, rel_min_height = 0.005,
                      color = alpha(1, .6)) +
  geom_vline(xintercept = 0, linetype = 'dashed', alpha = .1) +
  facet_wrap(~ Response, scales = 'free',
             shrink = TRUE, ncol = 1) +
  theme(axis.text.y = element_text(size = 7)) +
  scale_fill_gradient2() +
  theme(legend.position = 'none',
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(color = 'grey95')) +
  ylab('') +
  xlab('')
ggsave('fig/fire-effs.pdf', width = 6, height = 7)

rm(beta_df)
gc()


# Partial effect plots ----------------------------------------------------
which_var <- c('ctmx', 'crmin', 'cpr', 'cvs', 'cpr12', 'chd')
complete_name <- c('tmmx', 'rmin', 'pr', 'vs', 'prev_12mo_precip', 'housing_density')
assert_that(length(which_var) == length(complete_name))

partial_effs <- list()
n_iter <- length(post$lp__)
unique_ers <- unique(st_covs$NA_L3NAME)
pb <- txtProgressBar(max = length(unique_ers), style = 3)
for (v in seq_along(which_var)) {
  print(paste('Processing', which_var[v]))
  for (i in seq_along(unique_ers)) {
    setTxtProgressBar(pb, i)
    subdf <- st_covs %>%
      filter(NA_L3NAME == unique_ers[i]) %>%
      mutate(row_id = 1:n())
    X_sub <- X[st_covs$NA_L3NAME == unique_ers[i], ]
    cols <- grepl(which_var[v], colnames(X_sub))
    if (which_var[v] == 'cpr') {
      cols <- grepl(which_var[v], colnames(X_sub)) &
        !grepl('cpr12', colnames(X_sub))
    }

    effects <- array(dim = c(nrow(X_sub), 2, 3)) # 2 responses, 3: med, lo, hi
    for (j in 1:nrow(X_sub)) {  # month j
      for (k in 1:2) {          # response k
        vals <- X_sub[j, cols] %*% t(post$beta[, k, cols])
        effects[j, k, 1] <- quantile(vals, .025)
        effects[j, k, 2] <- median(vals)
        effects[j, k, 3] <- quantile(vals, .975)
      }
    }
    partial_effs[[length(partial_effs) + 1]] <- effects %>%
      reshape2::melt(varnames = c('row_id', 'response', 'quantity')) %>%
      as_tibble %>%
      mutate(response = ifelse(response == 1, 'negbinom', 'zeroinfl'),
             quantity = case_when(.$quantity == 1 ~ 'lo',
                                  .$quantity == 2 ~ 'med',
                                  .$quantity == 3 ~ 'hi'),
             var = which_var[v]) %>%
      spread(quantity, value) %>%
      left_join(select(subdf, row_id, NA_L3NAME, ym))
  }
}
close(pb)


name_df <- tibble(covariate = which_var,
                  untrans_cov = complete_name,
                  fancy_name = c('Air temperature (C)',
                                 'Relative humidity (%)',
                                 'Precipitation: same month (mm)',
                                 'Wind speed (m/s)',
                                 'Precipitation: 12 month (mm)',
                                 'Housing density (units per km^2)'))

effect_plot_df <- list()
for (i in seq_along(which_var)) {
  effect_plot_df[[i]] <- partial_effs %>%
    bind_rows %>%
    filter(var == which_var[i]) %>%
    left_join(st_covs) %>%
    mutate(covariate = which_var[i]) %>%
    left_join(name_df) %>%
    filter(response == 'negbinom') %>%
    mutate(covariate = complete_name[i]) %>%
    rename(covariate_value = !!complete_name[i]) %>%
    select(covariate, covariate_value, med, lo, hi,
           NA_L3NAME, ym, NA_L1NAME, NA_L2NAME,
           fancy_name)
}

annotation_df <- tibble(
  untrans_cov = c('tmmx', 'tmmx'),
  covariate_value = c(-7, 13),
  med = c(2, 3.5),
  label = c('Great Plains', 'Eastern Temperate Forests'),
  NA_L1NAME = c('GREAT PLAINS', 'EASTERN TEMPERATE FORESTS')
)

better_name <- distinct(ecoregion_df, NA_L1NAME) %>%
  mutate(l1_er = tools::toTitleCase(tolower(NA_L1NAME))) %>%
  mutate(ifelse(l1_er == 'Southern Semi-Arid Highlands',
                'Southern Semi-arid Highlands',
                l1_er))

p <- effect_plot_df %>%
  bind_rows %>%
  mutate(covariate_value = case_when(
    .$covariate == 'tmmx' ~ covariate_value - 273.15,
    .$covariate == 'rmin' ~ covariate_value / 100,
    TRUE ~ covariate_value)) %>%
  left_join(better_name) %>%
  ggplot(aes(covariate_value, y = med)) +
  theme_minimal() +
  geom_ribbon(aes(ymin = lo, ymax = hi, group = NA_L3NAME),
              color = NA, alpha = .04) +
  geom_line(alpha = .5, aes(group = NA_L3NAME, color = l1_er)) +
  ylab('Partial effect') +
  # geom_text(data = annotation_df %>% left_join(name_df),
  #           aes(label = label, color = NA_L1NAME), size = 3.5) +
  facet_wrap(~fancy_name, scales = 'free', ncol = 2) +
  scale_color_gdocs('') +
  xlab('') +
  theme(panel.grid.minor = element_blank(),
        legend.position = 'none')
p
ggsave('fig/count-partial-effs.png', width = 6, height = 4.5)


# Expected values vs. covariate values ------------------------------------
mu_df <- post$mu_full %>%
  reshape2::melt(varnames = c('iter', 'dim', 'id')) %>%
  as_tibble() %>%
  spread(dim, value) %>%
  mutate(expected_value = exp(`1`) * plogis(`2`)) %>%
  group_by(id) %>%
  summarize(med = median(expected_value),
            lo = quantile(expected_value, .05),
            hi = quantile(expected_value, .95))

st_covs <- st_covs %>%
  group_by(NA_L3NAME) %>%
  mutate(mean_hd = mean(housing_density),
         facet_factor = paste(NA_L2CODE, NA_L3NAME, sep = ': '),
         l2_er = tools::toTitleCase(tolower(as.character(NA_L2NAME))),
         l2_er = gsub(' and ', ' & ', l2_er),
         l2_er = gsub('Usa ', '', l2_er))
