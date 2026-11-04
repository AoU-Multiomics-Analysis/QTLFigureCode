# %% Locus and sequence-prediction methods
#| Define variant matching, local gene-track selection, PIP plotting, and PromoterAI window construction.
# @include notebooks/current/plots/Fig5_locus.R

# %% Read ACVRL1 eQTL and prediction result tables
#| Use the combined expression fine-mapping result, public cCRE table and GENCODE gene/exon rows, and the supplied PromoterAI predictions.
all_eqtl <- as_tibble(data.table::fread(susie_file("COMB.expression_SusieMerged.annotated.ancestry_skew.tsv.gz"),
  select = c("variant", "gene_name", "position", "pip", "cs_id", "gvs_max_af", "gvs_max_subpop", "log2_aFC"))) %>%
  mutate(rare = gvs_max_af < 0.01, chrom = str_extract(variant, "^chr[^_]+"))
eqtl <- all_eqtl %>% filter(gene_name == "ACVRL1")
qtl <- .acvrl1_qtl(eqtl)
frequencies <- build_variant_freq(all_eqtl) %>%
  semi_join(qtl %>% transmute(variant_key = normalize_variant_id(variant)), by = "variant_key")
promoter <- read_tsv(source_file("Fig5/reference/acvrl1_promoterAI.tsv"), show_col_types = FALSE)
ccre <- read_tsv(source_file("Fig5/reference/acvrl1_ccre_region.tsv"), show_col_types = FALSE)
# Calculate the same padded regional window used by the PIP plot.
regional <- regional_window_bp(qtl$position, qtl$pip >= 0.9, 51912316, 500) / 1e6
regional <- regional + c(-1, 1) * diff(regional) * 0.015
models <- s31_read_gene_models(input_file("gencode.v48.annotation.gtf.gz"), "chr12",
  regional[1] * 1e6, regional[2] * 1e6)
gene_spans <- models$spans
gene_exons <- models$exons

# %% Select the locus and build the regional PIP plot
#| Use PIP >=0.9 for the locus markers. Select the regional range, the promoter inset, overlapping cCREs, genes, and exons from their coordinates.
palette <- cs_palette(qtl$cs_id)
locus <- fig5_plot_acvrl1_locus(qtl, ccre)

# %% Join causal variants to the prediction windows
#| Match the transcript and allele identifiers. Keep the original 20-bp windows and signed-square prediction axis.
promoter_locus <- prepare_promoterAI_locus(promoter, qtl, "ACVRL1", "ENST00000419526.6",
  focal_position = 51912316, window_bp = 100)
facets <- build_causal_promoterAI_facets(promoter_locus, palette, flank_bp = 20, freq_lookup = frequencies)
predictions <- plot_causal_promoterAI_facets(promoter_locus, palette, flank_bp = 20,
  effect_threshold = 0.5, freq_lookup = frequencies,
  motif_annotations = tibble(facet_id = "chr12_51912316_C_T", xmin = -3.5, xmax = 5.5, tf = "SP1"),
  nucleotide_cols = nucleotide_colors)
stopifnot(nrow(predictions$data) == nrow(facets), n_distinct(facets$facet_id) == 3L)
# %% Draw fine-mapping and prediction panels
#| Panel A shows hg38 locus PIPs, credible sets, rare-variant markers, the promoter inset, ENCODE elements, and gene models. Panel B shows PromoterAI single-nucleotide predictions in three 20-bp windows. The signed-square prediction axis and focal-variant labels remain unchanged. The SP1 band is a manual annotation.
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
