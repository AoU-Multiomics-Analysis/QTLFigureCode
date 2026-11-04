# Plot builders use source-derived frames only. They do not read saved objects.
fig3_notebook_theme <- function() {
  theme_classic(base_family = "sans", base_size = 8) +
    theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
          strip.text = element_text(size = 8), strip.background = element_blank(),
          plot.title = element_blank(), plot.subtitle = element_blank(),
          plot.margin = margin(3, 3, 3, 3, "mm"))
}

fig3_panel_a <- function(layer_counts) {
  shared <- sum(layer_counts$prop[layer_counts$NumColocs > 0])
  layer_counts %>%
    mutate(NumColocs = factor(NumColocs, levels = 0:3)) %>%
    ggplot(aes(NumColocs, prop, fill = NumColocs)) +
    geom_col(width = .75) +
    scale_fill_manual(values = c("0" = "#D9D9D9", "1" = "#A9C1D9",
                                 "2" = "#7B9AB8", "3" = "#577590"), guide = "none") +
    annotate("segment", x = 1.5, xend = 4.5, y = .455, yend = .455,
             colour = "grey50", linewidth = .25) +
    annotate("text", x = 3, y = .495, size = 2.1, colour = "grey20",
             label = sprintf("%.1f%% shared across ≥1\nmolecular layer", 100 * shared)) +
    scale_y_continuous(labels = scales::label_percent(accuracy = 1), breaks = seq(0, .5, .1)) +
    coord_cartesian(ylim = c(0, .55)) +
    labs(x = "Number of colocalized\nmolecular layers", y = "Percent of genes") +
    fig3_notebook_theme() + theme(axis.line.x = element_blank(), axis.ticks.x = element_blank())
}

fig3_panel_b <- function(summaries) {
  main <- ggplot(summaries$signal_counts, aes(factor(n), nn / sum(nn), fill = factor(n))) +
    geom_col(width = .75) +
    scale_fill_manual(values = c("1" = "#D8BBBB", "2" = "#C89A9A", "3" = "#B37979",
                                 "4" = "#9E5E5E", "5" = "#8A4D4D", "6" = "#763F3F",
                                 "7" = "#643232"), guide = "none") +
    scale_y_continuous(labels = scales::label_percent(), expand = expansion(mult = c(0, .05))) +
    labs(x = "Number of\ncolocalized signals", y = "Percent of genes") +
    fig3_notebook_theme()
  inset <- ggplot(summaries$coloc_rate, aes(coloc_rate, reorder(type, coloc_rate), fill = type)) +
    geom_col(show.legend = FALSE) +
    scale_fill_manual(values = c(eQTL = "#007C8D", pQTL = "#D44258", sQTL = "#EFB144")) +
    scale_x_continuous(labels = scales::label_percent(accuracy = 1), breaks = c(0, .15, .3),
                       expand = expansion(mult = c(0, .05))) +
    labs(x = "% CS coloc", y = NULL) + fig3_notebook_theme() +
    theme(axis.text = element_text(size = 6), axis.title = element_text(size = 7),
          plot.margin = margin(0, 0, 0, 0))
  wrap_elements(full = main + inset_element(inset, left = .48, bottom = .60,
                                            right = .98, top = .98, ignore_tag = TRUE))
}

fig3_panel_c <- function(interaction_results) {
  interaction_results %>% filter(str_detect(term, ":")) %>%
    ggplot(aes(estimate, type)) +
    geom_linerange(aes(xmin = conf.low, xmax = conf.high), colour = "grey40") +
    geom_point(size = 2.8, colour = "#364B63") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    labs(x = "Effect of -log(MAF) on\ncross-ome correlation", y = NULL) +
    fig3_notebook_theme() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())
}

fig3_panel_d <- function(effects) {
  slopes <- effects$slopes %>%
    mutate(label = paste0("β = ", round(slope, 2), "\nN = ", scales::comma(n)))
  ggplot(effects$data, aes(posterior_mean.x, posterior_mean.y)) +
    geom_hline(yintercept = 0, linewidth = .25, colour = "grey85") +
    geom_vline(xintercept = 0, linewidth = .25, colour = "grey85") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = .35, colour = "grey65") +
    geom_point(alpha = .6, size = .9, colour = "#6EC5C3") +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = .7,
                colour = "red", fill = "red", alpha = .12) +
    geom_label(data = slopes, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
               hjust = -.05, vjust = 1.15, size = 2.4, linewidth = 0, fill = "white", alpha = .85) +
    facet_wrap(~maf_bin, nrow = 1) +
    coord_cartesian(xlim = c(-effects$limit, effects$limit), ylim = c(-effects$limit, effects$limit)) +
    labs(x = "Expression posterior", y = "Protein posterior") +
    fig3_notebook_theme() + theme(panel.spacing.x = unit(.8, "lines"))
}

fig3_panel_e <- function(enrichment, min_n = PANEL_E_MIN_N) {
  d <- as_tibble(enrichment) %>% filter(padj < .05, n >= min_n)
  unlabelled <- setdiff(d$feature, names(PANEL_E_LABELS))
  if (length(unlabelled)) message("Significant annotations without a display label: ",
                                  paste(unlabelled, collapse = ", "))
  d <- d %>% filter(feature %in% names(PANEL_E_LABELS)) %>%
    mutate(dir = if_else(fold > 1, "rare", "common")) %>% arrange(fold) %>%
    mutate(label = factor(PANEL_E_LABELS[feature], levels = PANEL_E_LABELS[feature]))
  if (!nrow(d)) return(ggplot() + theme_void())
  colours <- c(rare = "#8A4D4D", common = "#577590")
  limits <- c(min(d$lo) * .85, max(d$hi) * 1.25)
  breaks <- c(.5, .7, 1, 2, 5, 10, 20)
  breaks <- breaks[breaks > limits[1] & breaks < limits[2]]
  lx <- function(f) 10^(f * (log10(limits[2]) - log10(limits[1])))
  y_annotation <- nrow(d) + .85
  ggplot(d, aes(fold, label)) +
    annotate("rect", xmin = limits[1], xmax = 1, ymin = -Inf, ymax = Inf,
             fill = "grey95", colour = NA) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_linerange(aes(xmin = lo, xmax = hi), colour = "grey40") +
    geom_point(aes(colour = dir), size = 2.8) +
    annotate("segment", x = lx(-.004), xend = lx(-.055), y = y_annotation, yend = y_annotation,
             arrow = arrow(length = unit(5, "pt"), type = "closed"), colour = colours[["common"]], linewidth = .4) +
    annotate("segment", x = lx(.004), xend = lx(.055), y = y_annotation, yend = y_annotation,
             arrow = arrow(length = unit(5, "pt"), type = "closed"), colour = colours[["rare"]], linewidth = .4) +
    annotate("text", x = lx(-.065), y = y_annotation, label = "common", hjust = 1,
             size = 2.4, colour = colours[["common"]]) +
    annotate("text", x = lx(.065), y = y_annotation, label = "rare", hjust = 0,
             size = 2.4, colour = colours[["rare"]]) +
    scale_colour_manual(values = colours, guide = "none") +
    scale_x_log10(breaks = breaks, labels = breaks, limits = limits, expand = expansion(mult = 0)) +
    scale_y_discrete(expand = expansion(add = c(.6, 1.5))) +
    coord_cartesian(clip = "off") +
    labs(x = "MAF fold change vs.\nother colocalizing variants", y = NULL) +
    fig3_notebook_theme() +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          plot.margin = margin(l = 12, r = 24, t = 12, b = 6, unit = "pt"))
}

fig3_assemble <- function(panels) {
  coloc_counts <- panels$A + panels$B + plot_layout(widths = c(1, 1))
  top_row <- (coloc_counts | panels$C) + plot_layout(widths = c(2.5, 1))
  bottom_row <- (panels$D | panels$E) + plot_layout(widths = c(1.3, 1))
  (top_row / bottom_row) + plot_annotation(tag_levels = "A", theme = theme()) &
    theme(plot.tag = element_text(size = 11, face = "bold"))
}
