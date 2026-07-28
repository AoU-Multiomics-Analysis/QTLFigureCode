# =============================================================================
# Supplemental — full colocalization feature enrichment. NOT part of Figure 3.
# Produces FigS_ColocEnrichment.pdf:
#   analysis A — what makes a variant shared (coloc vs non-coloc, within pair)
#   analysis B — what makes a shared variant rare (within Coloc == TRUE)
#
# Panel E of the main figure shows only analysis A, and only the two
# well-powered pairs. Everything else lives here.
#
# Run from the project root:  Rscript fig3/supp_coloc_enrichment.R
# =============================================================================
source("fig3/00_setup.R")

enr <- coloc_feature_enrichment()

ANALYSIS_LABELS <- c(A = "A — shared vs layer-specific",
                     B = "B — rare vs common, among shared")

forest <- function(d, title) {
  ggplot(d, aes(x = odds_ratio, y = reorder(feature, odds_ratio),
                colour = flagged)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
    geom_pointrange(aes(xmin = conf.low, xmax = conf.high), size = 0.35) +
    facet_wrap(~pair, nrow = 1) +
    scale_x_continuous(trans = "log10") +
    # Grey = the pair was flagged underpowered for this analysis.
    scale_colour_manual(values = c(`FALSE` = "#364B63", `TRUE` = "grey70"),
                        guide = "none") +
    labs(x = "Adjusted odds ratio", y = NULL, title = title) +
    theme(plot.title = element_text(size = SZ_STRIP, face = "bold"))
}

figS <- forest(dplyr::filter(enr, analysis == "A"), ANALYSIS_LABELS[["A"]]) /
        forest(dplyr::filter(enr, analysis == "B"), ANALYSIS_LABELS[["B"]])

ggsave(file.path(OUTDIR, "FigS_ColocEnrichment.pdf"), figS,
       width = 14, height = 12, dpi = DPI, limitsize = FALSE)
ggsave(file.path(OUTDIR, "FigS_ColocEnrichment.png"), figS,
       width = 14, height = 12, dpi = DPI, limitsize = FALSE)
message("Assembled: ", file.path(OUTDIR, "FigS_ColocEnrichment.pdf"))
