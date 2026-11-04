# %% Read the saved plot inputs
#| Read the allele-frequency, ancestry-enrichment, allelic-effect, outlier, and gene-burden inputs.
# Reflow the figure around readable plot areas at the final print width.
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggbeeswarm)
  library(ggh4x)
  library(ggdist)
  library(grid)
})

source <- readRDS(paper_file('output/figure2_binned_draft/science_print_size/Science_display_plots.rds'))
showtext::showtext_auto(FALSE)
pt <- ggplot2::.pt
base_theme <- theme_classic(base_size = 8, base_family = 'Helvetica') +
  theme(text = element_text(family = 'Helvetica', size = 8),
        axis.title = element_text(size = 8),
        axis.title.x = element_text(margin = margin(t = 4)),
        axis.title.y = element_text(margin = margin(r = 4)),
        axis.text = element_text(size = 7, color = 'black'),
        axis.ticks = element_line(linewidth = .18),
        axis.ticks.length = unit(1, 'mm'),
        axis.line = element_line(linewidth = .18),
        strip.text = element_text(size = 8, margin = margin(3, 0, 3, 0)),
        strip.background = element_blank(),
        legend.text = element_text(size = 7), legend.title = element_blank(),
        panel.grid = element_blank(), plot.margin = margin(2, 2, 2, 2))
theme_set(base_theme)
plots <- map(source, ~.x + base_theme)
for (letter in names(plots)) {
  for (i in seq_along(plots[[letter]]$layers)) {
    if (inherits(plots[[letter]]$layers[[i]]$geom, 'GeomText')) {
      plots[[letter]]$layers[[i]]$aes_params$size <- 7 / pt
      plots[[letter]]$layers[[i]]$aes_params$family <- 'Helvetica'
    }
  }
}

# %% Four frequency bins
#| The unit is a high-PIP variant-trait pair. Keep <1%, 1-<5%, 5-<10%, and 10-50% bins and the recorded percentages.
plots$A <- plots$A +
  labs(x = 'Population-max MAF (%)', y = 'Causal variants (%)') +
  scale_x_discrete(labels = c('<1%', '1–<5%', '5–<10%', '10–50%'),
                    drop = FALSE, expand = expansion(add = .5)) +
  theme(axis.text.x = element_text(size = 7, angle = 0, hjust = .5),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = 'inside', legend.position.inside = c(.02, 1),
        legend.justification = c(0, 1), legend.direction = 'horizontal',
        legend.key.size = unit(3, 'mm'), legend.spacing.x = unit(1, 'mm'),
        legend.margin = margin(0, 0, 0, 0))

plots$B$layers[[2]]$aes_params$size <- .20
# Remove the common/rare annotation at the user's request.
plots$B$layers[[7]]$data <- plots$B$layers[[7]]$data %>% slice(0)
# %% Ancestry enrichment
#| Retain the saved enrichment tests and rare/common split. Display OTH as REM.
plots$B <- plots$B +
  labs(x = 'Ancestry group',
       y = 'Population enrichment') +
  scale_x_discrete(labels = function(x) stringr::str_replace(x, '^OTH$', 'REM')) +
  scale_y_continuous(transform = 'log10', breaks = 10^(0:4),
                     labels = expression(10^0, 10^1, 10^2, 10^3, 10^4)) +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_text(margin = margin(t = 4)),
        axis.text.y = element_text(size = 10))

plots$C$layers[[2]]$aes_params$point_size <- 1.4
plots$C$layers[[3]]$data <- plots$C$layers[[3]]$data %>%
  mutate(label = str_replace(label, fixed('+-'), '−'), x = 2.04)
# %% Allelic fold change
#| Keep the recorded distributions and test results. Zoom to 0-2 without removing source values.
plots$C <- plots$C +
  labs(x = '|log2 aFC|') +
  scale_x_continuous(breaks = seq(0, 2, .5), expand = expansion(mult = 0)) +
  coord_cartesian(xlim = c(0, 2), clip = 'on') +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = margin(2, 32, 2, 2))

plots$D$layers[[2]]$aes_params$size <- .10
# Keep annotation coordinates to preserve the axis ranges, but remove the text.
plots$D$layers[[7]]$data <- plots$D$layers[[7]]$data %>% mutate(label = '')
plots$D$facet$params$labeller <- as_labeller(c(
  'Expression\ndecrease' = 'Decreasing variants',
  'Expression\nincrease' = 'Increasing variants'))
# %% Cumulative effects
#| Keep the original expression cap below 400% and the pooled 4+ category. Thresholds are 60% and 150%.
plots$D <- plots$D +
  labs(x = 'Independent variants included',
       y = 'Expected heterozygous expression (%)') +
  theme(legend.position = 'none', axis.line.x = element_blank(),
        axis.ticks.x = element_blank(), panel.spacing.y = unit(3, 'mm'),
        plot.margin = margin(3, 2, 2, 2))

endpoints <- as_tibble(plots$E$data) %>%
  filter((type == 'Down' & focal_lower_bound == -100) |
           (type == 'Up' & focal_lower_bound == 175))
plots$E$layers[[4]]$aes_params$size <- 1.6
plots$E$layers[[5]]$data$y <- endpoints %>% filter(type == 'Down') %>% pull(odds_ratio) / 6
plots$E$layers[[6]]$data$y <- endpoints %>% filter(type == 'Up') %>% pull(odds_ratio) / 6
for (i in c(5, 6)) {
  plots$E$layers[[i]]$aes_params$hjust <- .5
  plots$E$layers[[i]]$aes_params$vjust <- .5
}
# F uses equal-width categories. Use those same positions for E's bins.
# E omits the central (-10, 10) reference bin, so that position stays empty.
bin_lower <- sort(unique(plots$F$data$lower))
bin_positions <- seq_along(bin_lower)
bin_transform <- local({
  boundaries <- bin_lower
  positions <- bin_positions
  scales::trans_new(
    'expression_bin_centers',
    transform = function(x) stats::approx(boundaries, positions, xout = x, rule = 2)$y,
    inverse = function(x) stats::approx(positions, boundaries, xout = x, rule = 2)$y,
    domain = range(boundaries))
})
plots$E <- plots$E +
  labs(y = 'Outlier\nenrichment') +
  scale_x_continuous(transform = bin_transform, limits = range(bin_lower),
                     expand = expansion(add = .6), breaks = NULL) +
  scale_y_continuous(transform = 'log2', breaks = 2^c(-4, 0, 4, 8),
                     labels = expression(2^"−4", 2^0, 2^4, 2^8)) +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        plot.margin = margin(3, 2, 1, 2))

# %% Gene burden
#| Retain aggregate gene-count medians and intervals. Participant rows are not printed.
plots$F <- plots$F +
  labs(x = 'Expression change bin (%)', y = 'Median genes\nper individual') +
  scale_x_discrete(labels = identity) +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 7, angle = 0, hjust = .5,
                                   lineheight = 1.1),
        strip.text = element_text(size = 7), plot.margin = margin(3, 2, 2, 2))
for (i in c(5, 6)) {
  plots$F$layers[[i]]$aes_params$y <- 8000
  plots$F$layers[[i]]$aes_params$vjust <- 1
}
for (i in c(8, 10)) {
  plots$F$layers[[i]]$aes_params$y <- 1700
}

checks <- imap_dfr(plots, function(p, letter) {
  stopifnot(identical(p$data, source[[letter]]$data))
  tibble(panel = letter, data_identical = TRUE,
         statistical_methods_identical = identical(map(p$layers, 'stat'),
                                                    map(source[[letter]]$layers, 'stat')))
})
stopifnot(all(checks$statistical_methods_identical))
write_tsv(checks, output_file('Fig2/data_preservation_checks.tsv'))

e_build <- ggplot_build(plots$E)
f_build <- ggplot_build(plots$F)
e_range <- e_build$layout$panel_params[[1]]$x.range
f_range <- f_build$layout$panel_params[[1]]$x.range
alignment <- tibble(lower_bound = sort(unique(plots$E$data$focal_lower_bound))) %>%
  mutate(bin_position = match(lower_bound, bin_lower),
         E_fraction = (bin_transform$transform(lower_bound) - e_range[[1]]) / diff(e_range),
         F_fraction = (bin_position - f_range[[1]]) / diff(f_range),
         aligned = abs(E_fraction - F_fraction) < 1e-10)
stopifnot(all(alignment$aligned))
write_tsv(alignment, output_file('Fig2/EF_alignment_check.tsv'))


# %% Render panel A
save_plot(plots$A, "Fig2/panel_A", 150, 105)
print(plots$A)

# %% Render panel B
save_plot(plots$B, "Fig2/panel_B", 150, 105)
print(plots$B)

# %% Render panel C
#| Clip the distributions at 2. Keep the statistical labels in the right
#| margin, as in the current vector export. Do not clip the entire panel grob.
render_panel_C <- function(plot) {
  result <- ggplot2::ggplotGrob(plot)
  index <- which(result$layout$name == 'panel')
  panel <- result$grobs[[index]]
  is_text <- vapply(panel$children, inherits, logical(1), what = 'text')
  stopifnot(sum(is_text) == 1)
  geometry <- do.call(grid::grobTree,
                     c(unname(as.list(panel$children[!is_text])),
                       list(vp = grid::viewport(clip = 'on'))))
  panel <- grid::setChildren(panel, do.call(grid::gList,
                            c(list(geometry), unname(as.list(panel$children[is_text])))))
  if (!is.null(panel$vp)) panel$vp$clip <- FALSE
  result$grobs[[index]] <- panel
  result$layout$clip[[index]] <- 'off'
  result
}
c_artwork <- render_panel_C(plots$C)
save_plot(c_artwork, "Fig2/panel_C", 150, 105)
grid::grid.newpage()
grid::grid.draw(c_artwork)

# %% Render panel D
save_plot(plots$D, "Fig2/panel_D", 150, 105)
print(plots$D)

# %% Outlier enrichment and burden
#| Align E and F on the same expression-bin centers. The central reference bin is omitted from E.
ef <- (plots$E / plots$F) + plot_layout(heights = c(24, 33))
save_plot(ef, "Fig2/panels_EF", 183, 90)
print(ef)

# %% Rebuild the compact layout
layout_plots <- plots
layout_plots$C <- wrap_elements(full = c_artwork)
fig2 <- wrap_plots(layout_plots, design = "AB\nCD\nEE\nFF", widths = c(81, 98), heights = c(48, 62, 24, 33)) + plot_annotation(tag_levels = "A")
save_plot(fig2, "Fig2/Figure2_rebuilt", 183, 188)
print(fig2)

# %% Reference figure
#| Complete figure layout.
show_reference(2)
