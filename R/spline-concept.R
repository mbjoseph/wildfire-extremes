library(patchwork)
library(splines)
library(tidyverse)

theme_set(theme_minimal() + 
            theme(panel.grid.minor = element_blank(), 
                  panel.grid.major = element_blank(), 
                  axis.text.x = element_blank(), 
                  plot.title = element_text(size=10)))

m <- 5
pal <- 2
axis_alpha <- .2

# Simulating function values ----------------------------------------------
n <- 100
x <- seq(0, 1, length.out = n)
unweighted_basis <- bs(x, df = 5, intercept = TRUE) %>%
  as.data.frame() %>%
  as_tibble %>%
  mutate(x = x) %>%
  gather(key = 'dim', 'basis_value', -x) %>%
  mutate(dim = as.integer(dim), 
         weight = 1,
         type = 'Unweighted')


unweight_plot <- unweighted_basis %>%
  ggplot(aes(x, basis_value, color = factor(dim))) + 
  geom_line() + 
  scale_color_brewer(palette = pal, type = 'qual') + 
  xlab('') + 
  ylab('') + 
  theme(legend.position = 'none') + 
  ylim(-1, 1) + 
  geom_vline(xintercept = 0, alpha = axis_alpha) + 
  geom_hline(yintercept = 0, alpha = axis_alpha) + 
  ggtitle('B-splines')
unweight_plot

# Weighted basis vectors --------------------------------------------------

weights <- c(-1, 1, -.5, 1, 1)

weighted_basis <- unweighted_basis %>%
  mutate(type = 'Weighted', 
         weight = weights[dim], 
         weighted_value = weight * basis_value)

function_val <- weighted_basis %>%
  group_by(x) %>%
  summarize(f_x = sum(weighted_value))

annotation_df <- weighted_basis %>%
  group_by(dim) %>%
  mutate(abs_weighted = abs(weighted_value), 
         is_max = abs_weighted == max(abs_weighted)) %>%
  filter(is_max) %>%
  ungroup

weight_plot <- weighted_basis %>%
  ggplot(aes(x)) + 
  geom_line(aes(y = weighted_value, color = factor(dim))) + 
  scale_color_brewer(palette = pal, type = 'qual') + 
  xlab('') +
  ylab('') + 
  theme(legend.position = 'none', 
        axis.text.y = element_blank()) + 
  geom_vline(xintercept = 0, alpha = axis_alpha) + 
  geom_hline(yintercept = 0, alpha = axis_alpha) + 
  ggtitle('Weighted B-splines') + 
  geom_text(aes(x = .04 + x * .93, 
                y = .55 * weighted_value, 
                label = weight, 
                color = factor(dim)), 
            data = annotation_df, size = 3) + 
  ylim(-1, 1)

function_plot <- function_val %>%
  ggplot(aes(x, f_x)) +
  xlab('') + 
  ylab('') + 
  geom_line() +  
  geom_vline(xintercept = 0, alpha = axis_alpha) + 
  geom_hline(yintercept = 0, alpha = axis_alpha) + 
  ggtitle('Function value') + 
  theme(axis.text.y = element_blank()) + 
  ylim(-1, 1)

p <- unweight_plot + weight_plot + function_plot
p

ggsave(filename = 'fig/figure_3.pdf', plot = p, width = 6, height = 2)
