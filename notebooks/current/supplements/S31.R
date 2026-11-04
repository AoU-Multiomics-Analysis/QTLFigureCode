# %% Read the saved ACVRL1 panels
#| Restore S31 from the saved Figure 5 ACVRL1 locus and PromoterAI objects. No fine-mapping, motif scan, or sequence model is run. The saved objects include the annotations and prediction coordinates needed to draw the panels.
locus <- readRDS(source_file("Fig5/objects/D.rds"))
predictions <- readRDS(source_file("Fig5/objects/E.rds"))
stopifnot(inherits(locus, "patchwork"), inherits(predictions, "ggplot"),
          all(c("promoterAI", "relative_position", "facet_id") %in% names(predictions$data)))

# %% Restore fine-mapping and prediction panels
#| Panel A shows saved hg38 locus PIPs, credible sets, rare-variant markers, the promoter inset, ENCODE elements, and gene models. Panel B shows saved PromoterAI single-nucleotide predictions in three 20-bp windows. The signed-square prediction axis and focal-variant labels remain unchanged. The SP1 band is a saved manual annotation.
remove_headings <- function(plot) {
  if (inherits(plot, "patchwork")) {
    plot$patches$annotation$title <- NULL
    plot$patches$annotation$subtitle <- NULL
    plot$patches$annotation$theme <- theme()
    plot$patches$plots <- lapply(plot$patches$plots, remove_headings)
  }
  plot + labs(title = NULL, subtitle = NULL) +
    theme(plot.title = element_blank(), plot.subtitle = element_blank())
}
locus <- remove_headings(locus)
# Use ASCII log2 in strip text. Some device fonts do not contain subscript 2.
# This changes the display labels only; all prediction values are retained.
original_labeller <- predictions$facet$params$labeller
predictions$facet$params$labeller <- function(labels, multi_line = TRUE) {
  labels <- lapply(labels, function(x) stringr::str_replace_all(as.character(x), fixed("\u2082"), "2"))
  original_labeller(labels, multi_line = multi_line)
}
class(predictions$facet$params$labeller) <- c("function", "labeller")
locus_panel <- patchwork::wrap_elements(full = patchwork::patchworkGrob(locus))
figure <- locus_panel / predictions + patchwork::plot_layout(heights = c(76, 138)) +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(plot.title = element_blank(), plot.subtitle = element_blank())
print(figure)

# %% Save aggregate outputs
#| Export the two-panel figure and counts per prediction window. Predictions describe model output, not an experimental test of regulatory activity or pathogenicity. No variant or participant source table is exported.
save_plot(figure, "supplements/S31", 184, 220)
readr::write_tsv(predictions$data %>% group_by(facet_id) %>%
                   summarise(substitution_predictions = n(), positions = n_distinct(relative_position),
                             minimum_prediction = min(promoterAI), maximum_prediction = max(promoterAI), .groups = "drop"),
                 output_file("supplements/S31_prediction_summary.tsv"))
