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

# enrichment_forest_panel() (Fig3/00_setup.R) renders a placeholder instead of
# faceting on an empty frame when a given analysis has zero rows -- see its
# doc comment for why that half can legitimately be empty even when the other
# is not.
figS <- enrichment_forest_panel(dplyr::filter(enr, analysis == "A"), ANALYSIS_LABELS[["A"]]) /
        enrichment_forest_panel(dplyr::filter(enr, analysis == "B"), ANALYSIS_LABELS[["B"]])

ggsave(file.path(OUTDIR, "FigS_ColocEnrichment.pdf"), figS,
       width = 14, height = 12, dpi = DPI, limitsize = FALSE)
ggsave(file.path(OUTDIR, "FigS_ColocEnrichment.png"), figS,
       width = 14, height = 12, dpi = DPI, limitsize = FALSE)
message("Assembled: ", file.path(OUTDIR, "FigS_ColocEnrichment.pdf"))
