source('R/thresholds/explore.R')
library(ismev)

# Evaluating stability of estimates to select a threshold ------------
n <- 500
threshold <- seq(mean(d$fire_size), quantile(d$fire_size, .99),
                  length.out = n)
params <- matrix(nrow = n, ncol = 2)
se <- matrix(nrow = n, ncol = 2)
sample_size <- rep(NA, n)
fits <- vector(mode = 'list', length = n)

range_df <- data.frame(threshold)

for (i in seq_along(threshold)) {
  fits[[i]] <- gpd.fit(d$fire_size, threshold = threshold[i])
  range_df$sigma[i] <- fits[[i]]$mle[1]
  range_df$ksi[i] <- fits[[i]]$mle[2]
  range_df$sigma_se[i] <- fits[[i]]$se[1]
  range_df$ksi_se[i] <- fits[[i]]$se[2]
  range_df$sample_size[i] <- fits[[i]]$nexc
}

range_df <- range_df %>%
  mutate(sigma_star = sigma - ksi * threshold)


# Graphically evaluate stability of estimates
# for a valid threshold, estimates of sigma_star and ksi should be stable
# above the threshold
ggplot(range_df, aes(x = threshold, y = ksi)) +
  geom_point(alpha = .5) +
  geom_segment(aes(xend = threshold,
                   y = ksi - ksi_se,
                   yend = ksi + ksi_se)) +
  theme_minimal() +
  ylab(expression(hat(xi)))

ggplot(range_df, aes(x = threshold, y = sigma_star)) +
  geom_point(alpha = .5) +
  theme_minimal() +
  ylab(expression(sigma^'*'))
