# =============================================================================
# Supplemental — cross-QTL colocalization QC. NOT part of Figure 3.
# Produces FigS_CrossQTLColocQC.pdf:
#   row 1: posterior distributions (GLCP / RCP / SCP) by assay pair
#   row 2: number of colocalizations discovered across posterior thresholds
#
# The original wrote this figure to CrossQTLPosteriorDistributions.pdf twice,
# with different content; the first save was silently overwritten.
#
# Run from the project root:  Rscript fig3/supp_coloc_qc.R
# =============================================================================
source("fig3/00_setup.R")

gene_lvl    <- coloc_gene_level()
cs_lvl      <- coloc_cs_level()
variant_lvl <- coloc_variant_level()

density_panel <- function(d, col, show_legend) {
  ggplot(d, aes(x = .data[[col]], fill = type)) +
    # The original applied alpha = 0.7 only to the FIRST density panel
    # (GLCP); RCP and SCP used geom_density()'s opaque default. Applying a
    # uniform alpha to all three here is a deliberate, intentional deviation
    # from the original -- uniform alpha reads better across the row -- not
    # an oversight.
    geom_density(alpha = 0.7, show.legend = show_legend) +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.y = element_blank())
}

# wrap_elements() folds the three density panels into a single taggable
# element, so plot_annotation(tag_levels = "A") below sees this row as ONE
# leaf (tag "A") rather than walking into it and tagging each density panel
# separately (which would shift the discovery plot's tag from "B" to "D").
posterior_dists <- wrap_elements(full = density_panel(gene_lvl,    "GLCP", FALSE) +
                                         density_panel(cs_lvl,      "RCP",  FALSE) +
                                         density_panel(variant_lvl, "SCP",  TRUE))

THRESHOLDS <- c(0.01, 0.1, 0.5, 0.9)

discovery <- function(d, col, label) {
  purrr::map_dfr(THRESHOLDS, function(t) {
    d %>%
      group_by(type) %>%
      summarise(n_tested = dplyr::n(),
                value    = sum(.data[[col]] > t, na.rm = TRUE),
                .groups  = "drop") %>%
      mutate(threshold = t, posterior = label)
  })
}

discovery_by_type <- bind_rows(
  discovery(gene_lvl,    "GLCP", "GLCP"),
  discovery(cs_lvl,      "RCP",  "RCP"),
  discovery(variant_lvl, "SCP",  "SCP")
)

discovery_plot <- discovery_by_type %>%
  ggplot(aes(x = threshold, y = log10(value), group = type, color = type)) +
    geom_point() +
    geom_line() +
    facet_wrap(~posterior, nrow = 1) +
    labs(x = "Posterior threshold", y = "log10(Number Colocs)")

SuppColocQC <- posterior_dists / discovery_plot +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = SZ_TAG, face = "bold"))

ggsave(file.path(OUTDIR, "FigS_CrossQTLColocQC.pdf"), SuppColocQC,
       height = 7, width = 10, dpi = DPI)
message("Wrote ", file.path(OUTDIR, "FigS_CrossQTLColocQC.pdf"))
