discovery_plot <- function(panel, y_label) {
  ggplot(panel$data, aes(sample_size, credible_sets)) +
    geom_polygon(data = panel$envelope$boundary, aes(sample_size, credible_sets),
      inherit.aes = FALSE, fill = "#E4E6E8", colour = "#B9BEC2", linewidth = 0.35) +
    geom_line(data = panel$curve, colour = "#667780", linewidth = 0.45) +
    geom_point(aes(fill = colour, shape = marker), size = 2.3, colour = "#263B45", stroke = 0.3) +
    ggrepel::geom_text_repel(aes(label = label), size = 2.4, colour = "#263B45",
      segment.colour = "#87949A", segment.size = 0.25, seed = 14, max.overlaps = Inf) +
    annotate("segment", x = max(panel$gtex$sample_size), xend = max(panel$gtex$sample_size) * 1.9,
      y = median(panel$gtex$credible_sets), yend = median(panel$gtex$credible_sets),
      colour = "#87949A", linewidth = 0.3) +
    annotate("text", x = max(panel$gtex$sample_size) * 2.0, y = median(panel$gtex$credible_sets),
      label = "GTEx v10", size = 2.4, hjust = 0, colour = "#66737A") +
    scale_shape_identity() + scale_fill_identity() +
    scale_x_continuous(labels = scales::label_comma(), expand = expansion(mult = c(.04, .10))) +
    scale_y_continuous(labels = scales::label_comma(), expand = expansion(mult = c(.04, .12))) +
    labs(x = "eQTL sample size", y = y_label) +
    theme_classic(base_size = 8) +
    theme(legend.position = "none", axis.text = element_text(colour = "#263B45"),
      axis.title = element_text(colour = "#263B45"), plot.margin = margin(5, 8, 5, 5))
}
