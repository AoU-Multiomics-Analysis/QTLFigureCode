# %% Figure 2 plot definitions
#| These functions construct the plots from the results calculated in this notebook.
#| No saved plot supplies data, models, or plot layers.
fig2_theme <- function() {
  theme_classic(base_size = 8, base_family = "Helvetica") +
    theme(text = element_text(family = "Helvetica", size = 8),
          axis.title = element_text(size = 8),
          axis.title.x = element_text(margin = margin(t = 4)),
          axis.title.y = element_text(margin = margin(r = 4)),
          axis.text = element_text(size = 7, color = "black"),
          axis.ticks = element_line(linewidth = .18),
          axis.ticks.length = grid::unit(1, "mm"),
          axis.line = element_line(linewidth = .18),
          strip.text = element_text(size = 8, margin = margin(3, 0, 3, 0)),
          strip.background = element_blank(),
          legend.text = element_text(size = 7), legend.title = element_blank(),
          panel.grid = element_blank(), plot.margin = margin(2, 2, 2, 2))
}

# Keep the current Figure 1 ancestry colors. ADMIX and OTH are one source group.
fig2_ancestry_colors <- function() {
  palette <- rcartocolor::carto_pal(8, "Temps")
  c(OTH = palette[1], ADMIX = palette[1], EUR = palette[2], AMR = palette[3],
    AFR = palette[4], SAS = palette[5], EAS = palette[6], MID = palette[7], AoU = palette[8])
}

fig2_plot_a <- function(counts) {
  bins <- c("<1%", "1–<5%", "5–<10%", "10–50%")
  counts <- counts %>% mutate(Assay = factor(Assay, c("Expression", "Proteomics", "Splicing")),
                              frequency_bin = factor(frequency_bin, bins))
  bin_lines <- tibble(x = seq_along(bins) - .4, xend = seq_along(bins) + .4)
  ggplot(counts, aes(frequency_bin, percent, color = Assay, shape = Assay)) +
    geom_vline(xintercept = 1.5, linetype = "dashed", color = "grey50", linewidth = .18) +
    annotate("text", x = 1.58, y = 63, label = "GTEx", hjust = 0,
             family = "Helvetica", size = 7 / ggplot2::.pt) +
    annotate("segment", x = 1.58, xend = 2.03, y = 59.5, yend = 59.5,
             linewidth = .18, arrow = grid::arrow(length = grid::unit(1.2, "mm"), type = "closed")) +
    geom_segment(data = bin_lines, aes(x, xend = xend, y = 0, yend = 0),
                 inherit.aes = FALSE, color = "grey60", linewidth = .18) +
    geom_point(position = position_dodge(width = .9), size = 2.1, stroke = .2) +
    scale_x_discrete(drop = FALSE, expand = expansion(add = .5)) +
    scale_color_manual(values = c(Expression = "#3097A6", Proteomics = "#DF687B", Splicing = "#DFA545"),
                       labels = c(Expression = "eQTL", Proteomics = "pQTL", Splicing = "sQTL")) +
    scale_shape_manual(values = c(Expression = 16, Proteomics = 18, Splicing = 15),
                       labels = c(Expression = "eQTL", Proteomics = "pQTL", Splicing = "sQTL")) +
    scale_y_continuous(limits = c(0, 75), breaks = seq(0, 60, 20),
                       labels = scales::label_number(suffix = "%"),
                       expand = expansion(mult = c(.03, .01))) +
    labs(x = "Population-max MAF (%)", y = "Causal variants (%)") + fig2_theme() +
    theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "inside", legend.position.inside = c(.02, 1),
          legend.justification = c(0, 1), legend.direction = "horizontal",
          legend.key.size = grid::unit(3, "mm"), legend.spacing.x = grid::unit(1, "mm"),
          legend.margin = margin(0, 0, 0, 0))
}

fig2_plot_b <- function(ancestry) {
  pd <- position_dodge(width = .8)
  ggplot(ancestry$skew_df, aes(pop, odds_ratio, fill = pop, color = pop,
                              group = interaction(pop, rare_lab))) +
    geom_hline(yintercept = 5, linetype = "dashed", linewidth = .18, color = "grey40") +
    ggbeeswarm::geom_quasirandom(width = .15, dodge.width = .8, size = .20, stroke = 0,
                                alpha = .7, show.legend = FALSE) +
    geom_violin(width = .9, alpha = .3, position = pd, color = NA, trim = FALSE, show.legend = FALSE) +
    stat_summary(fun = median, geom = "crossbar", linewidth = .18, width = .1,
                 alpha = .6, position = pd, color = "black", show.legend = FALSE) +
    scale_color_manual(values = fig2_ancestry_colors()) +
    scale_fill_manual(values = fig2_ancestry_colors()) +
    geom_text(data = ancestry$skew_counts, aes(pop, 9000, label = scales::comma(n)),
              inherit.aes = FALSE, size = 7 / ggplot2::.pt, family = "Helvetica",
              color = "grey40", fontface = "bold") +
    geom_segment(data = ancestry$pair_line_df, aes(x = x, xend = xend, y = y, yend = yend),
                 inherit.aes = FALSE, color = "grey60", linewidth = .18) +
    annotate("segment", x = .85, xend = 7.15, y = 18000, yend = 18000,
             linewidth = .18, color = "grey40") +
    annotate("text", x = 4, y = 26000, label = "Ancestry-enriched variants",
             size = 7 / ggplot2::.pt, family = "Helvetica", fontface = "bold", color = "grey40") +
    scale_x_discrete(labels = function(x) str_replace(x, "^OTH$", "REM")) +
    scale_y_continuous(transform = "log10", breaks = 10^(0:4),
                       labels = expression(10^0, 10^1, 10^2, 10^3, 10^4)) +
    labs(x = "Ancestry group", y = "Population enrichment") +
    coord_cartesian(clip = "off") + fig2_theme() +
    theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.x = element_text(margin = margin(t = 4)), axis.text.y = element_text(size = 10))
}

fig2_plot_c <- function(afc) {
  regression <- afc$log2_afc_regression %>% mutate(label = str_replace(label, fixed("+-"), "−"), x = 2.04)
  afc$afc_plot_df %>% filter(!str_detect(color_group, "NA")) %>%
    ggplot(aes(abs_log2_aFC, name, color = color_group, fill = color_group)) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = .18) +
    ggdist::stat_halfeye(fill_type = "segments", point_interval = "median_qi", show.legend = FALSE,
                         height = 1, point_size = 1.4, interval_size = .18,
                         slab_linewidth = .18, slab_alpha = .7) +
    scale_fill_manual(values = afc$afc_palette) + scale_color_manual(values = afc$afc_palette) +
    geom_text(data = regression, aes(x, name, label = label), inherit.aes = FALSE,
              hjust = 0, size = 7 / ggplot2::.pt, family = "Helvetica") +
    labs(x = "|log2 aFC|", y = NULL) +
    scale_x_continuous(breaks = seq(0, 2, .5), expand = expansion(mult = 0)) +
    coord_cartesian(xlim = c(0, 2), clip = "on") + fig2_theme() +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          plot.margin = margin(2, 32, 2, 2))
}

fig2_plot_d <- function(display) {
  band_df <- tibble(direction = c("Expression\ndecrease", "Expression\nincrease"),
                    xmin = -Inf, xmax = Inf, ymin = c(.45, 1.45), ymax = c(.6, 4))
  # Empty labels retain the current annotation coordinates and axis extents.
  band_label_df <- tibble(direction = c("Expression\ndecrease", "Expression\nincrease"),
                          x = c(.5, 4.45), y = c(.43, 4.15), label = "", hjust = c(0, 1))
  labels <- display$label_df %>% mutate(label = scales::percent(prop, accuracy = 1),
                                       y = if_else(direction == "Expression\ndecrease", 1.1, 4.65))
  ggplot(display$plot_df, aes(factor(k_label), expected_het_expr, color = direction, fill = direction)) +
    geom_violin(width = .8, alpha = .3, color = NA, trim = FALSE, show.legend = FALSE) +
    ggbeeswarm::geom_quasirandom(width = .1, size = .10, stroke = 0, alpha = .2, show.legend = FALSE) +
    geom_hline(yintercept = 1, linetype = 2, color = "grey40", linewidth = .18) +
    facet_wrap(~direction, scales = "free_y", nrow = 2, strip.position = "top",
               labeller = as_labeller(c("Expression\ndecrease" = "Decreasing variants",
                                        "Expression\nincrease" = "Increasing variants"))) +
    scale_fill_manual(values = c("#B65C73", "#2A9D8F")) +
    scale_color_manual(values = c("#B65C73", "#2A9D8F")) +
    stat_summary(fun = median, geom = "crossbar", linewidth = .18, width = .2,
                 alpha = .2, color = "black", show.legend = FALSE) +
    ggh4x::facetted_pos_scales(y = list(
      scale_y_continuous(breaks = c(.5, .7, .9), labels = scales::percent_format(accuracy = 1)),
      scale_y_continuous(breaks = c(1.5, 2.5, 3.5), labels = scales::percent_format(accuracy = 1)))) +
    geom_text(data = labels, aes(k_label, y, label = label), inherit.aes = FALSE,
              size = 7 / ggplot2::.pt, family = "Helvetica") +
    geom_rect(data = band_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE, fill = "grey40", alpha = .05) +
    geom_text(data = band_label_df, aes(x, y, label = label, hjust = hjust), inherit.aes = FALSE,
              size = 7 / ggplot2::.pt, family = "Helvetica", color = "grey35", fontface = "bold") +
    labs(x = "Independent variants included", y = "Expected heterozygous expression (%)") +
    coord_cartesian(clip = "off") + fig2_theme() +
    theme(legend.position = "none", axis.line.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title.y.left = element_text(margin = margin(r = 3)),
          panel.spacing.y = grid::unit(3, "mm"), plot.margin = margin(3, 2, 2, 2))
}

fig2_bin_transform <- function(lower) {
  boundaries <- sort(unique(lower))
  positions <- seq_along(boundaries)
  scales::trans_new("expression_bin_centers",
                    transform = function(x) stats::approx(boundaries, positions, xout = x, rule = 2)$y,
                    inverse = function(x) stats::approx(positions, boundaries, xout = x, rule = 2)$y,
                    domain = range(boundaries))
}

fig2_plot_e <- function(outliers, median) {
  data <- outliers %>% filter(focal_lower_bound != -10)
  endpoints <- data %>% filter((type == "Down" & focal_lower_bound == -100) |
                                (type == "Up" & focal_lower_bound == 175))
  outlier_pal <- c("#D44D5CFF", "#046E8FFF")
  ggplot(data, aes(focal_lower_bound, odds_ratio, color = type, fill = type, group = type)) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = .18) +
    geom_line(show.legend = FALSE, linewidth = .18) +
    geom_linerange(aes(ymin = ci_low, ymax = ci_high), show.legend = FALSE, linewidth = .18) +
    geom_point(size = 1.6, shape = 21, alpha = .7, stroke = .2, show.legend = FALSE) +
    scale_fill_manual(values = outlier_pal) + scale_color_manual(values = outlier_pal) +
    annotate("text", x = -100, y = endpoints$odds_ratio[endpoints$type == "Down"] / 6,
             label = "Underlier", color = "#D44D5CFF", size = 7 / ggplot2::.pt,
             family = "Helvetica", hjust = .5, vjust = .5) +
    annotate("text", x = 175, y = endpoints$odds_ratio[endpoints$type == "Up"] / 6,
             label = "Overlier", color = "#046E8FFF", size = 7 / ggplot2::.pt,
             family = "Helvetica", hjust = .5, vjust = .5) +
    scale_x_continuous(transform = fig2_bin_transform(median$lower), limits = range(median$lower),
                       expand = expansion(add = .6), breaks = NULL) +
    scale_y_continuous(transform = "log2", breaks = 2^c(-4, 0, 4, 8),
                       labels = expression(2^"−4", 2^0, 2^4, 2^8)) +
    labs(y = "Outlier\nenrichment") + coord_cartesian(clip = "off") + fig2_theme() +
    theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.x = element_blank(), axis.title.x = element_blank(),
          axis.text.y = element_text(size = 10), plot.margin = margin(3, 2, 1, 2))
}

fig2_plot_f <- function(median) {
  max_median <- max(median$q75_genes)
  one <- median %>% group_by(gene_type) %>% slice(1) %>% ungroup()
  # These annotations sum per-bin medians. They are not participant totals.
  tail_medians <- median %>% group_by(gene_type) %>%
    summarise(down = sum(median_genes[lower < -50]), up = sum(median_genes[lower >= 50]), .groups = "drop")
  ggplot(median, aes(reorder(pct_bin, lower), median_genes)) +
    annotate("rect", xmin = .5, xmax = 2.45, ymin = -1, ymax = max_median + 100, fill = "grey40", alpha = .05) +
    annotate("rect", xmin = 7.5, xmax = 13.45, ymin = -1, ymax = max_median + 100, fill = "grey40", alpha = .05) +
    geom_col(width = .8, alpha = .6, linewidth = .18, color = "#462255FF", fill = "#462255FF") +
    geom_errorbar(aes(ymin = q25_genes, ymax = q75_genes), alpha = .7, width = .2, linewidth = .18) +
    facet_wrap(~gene_type, nrow = 2, strip.position = "right") +
    geom_text(data = one, x = 1.5, y = 8000, label = "Del-Like", size = 7 / ggplot2::.pt,
              family = "Helvetica", color = "grey35", fontface = "bold", vjust = 1) +
    geom_text(data = one, x = 10.475, y = 8000, label = "Dup-Like", size = 7 / ggplot2::.pt,
              family = "Helvetica", color = "grey35", fontface = "bold", vjust = 1) +
    annotate("segment", x = .6, xend = 2.25, y = 1000, yend = 1000, color = "black", linewidth = .18) +
    geom_text(data = tail_medians, aes(label = paste(scales::comma(down), "Genes")),
              inherit.aes = FALSE, x = 1.425, y = 1700, size = 7 / ggplot2::.pt,
              family = "Helvetica", color = "grey35", fontface = "bold") +
    annotate("segment", x = 8, xend = 12.95, y = 1000, yend = 1000, color = "black", linewidth = .18) +
    geom_text(data = tail_medians, aes(label = paste(scales::comma(up), "Genes")),
              inherit.aes = FALSE, x = 10.475, y = 1700, size = 7 / ggplot2::.pt,
              family = "Helvetica", color = "grey35", fontface = "bold") +
    labs(x = "Expression change bin (%)", y = "Median genes\nper individual") + fig2_theme() +
    theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.x = element_text(size = 7, angle = 0, hjust = .5, lineheight = 1.1),
          strip.text = element_text(size = 7), plot.margin = margin(3, 2, 2, 2))
}

# %% Ancestry-frequency plot definition
#| OTH is displayed as REM. The source ancestry keys are retained.
fig2_plot_s16 <- function(bins) {
  stack <- bins$stack_df %>% mutate(pop = toupper(gvs_max_subpop),
                                    pop = if_else(pop == "ADMIX", "OTH", pop))
  ggplot(stack, aes(x, prop, fill = pop)) +
    geom_area(position = "stack", alpha = .95, linewidth = .5, na.rm = TRUE) +
    geom_smooth(inherit.aes = FALSE, data = bins$line_df, aes(x, prop_shared_smooth),
                color = "black", linewidth = .4, se = FALSE, na.rm = TRUE) +
    geom_vline(xintercept = .5, color = "white", linewidth = .5) +
    scale_fill_manual(values = fig2_ancestry_colors(), labels = function(x) if_else(x == "OTH", "REM", x)) +
    scale_x_log10(breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 5e-1),
                  labels = c("0.0001", "0.001", "0.01", "0.1", "0.5"), expand = c(0, 0)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       expand = c(0, .001), breaks = c(.25, .5, .75, 1)) +
    labs(x = "AoU allele frequency", y = "Proportion of ancestry-enriched variants", fill = "Enriched\nancestry") +
    fig2_theme() + theme(axis.line = element_line(color = "black", linewidth = .3),
                         axis.ticks = element_line(linewidth = .2),
                         axis.title.y = element_text(margin = margin(r = 20)),
                         legend.title = element_text(size = 8), legend.position = "inside",
                         legend.position.inside = c(.9, .5), plot.margin = margin(14, 32, 14, 14))
}

# %% aFC rendering definition
#| Clip the distributions at 2. Keep statistical labels in the right margin.
render_panel_C <- function(plot) {
  result <- ggplot2::ggplotGrob(plot)
  index <- which(result$layout$name == "panel")
  panel <- result$grobs[[index]]
  is_text <- vapply(panel$children, inherits, logical(1), what = "text")
  stopifnot(sum(is_text) == 1)
  geometry <- do.call(grid::grobTree, c(unname(as.list(panel$children[!is_text])),
                                      list(vp = grid::viewport(clip = "on"))))
  panel <- grid::setChildren(panel, do.call(grid::gList,
                                           c(list(geometry), unname(as.list(panel$children[is_text])))))
  if (!is.null(panel$vp)) panel$vp$clip <- FALSE
  result$grobs[[index]] <- panel
  result$layout$clip[[index]] <- "off"
  result
}
