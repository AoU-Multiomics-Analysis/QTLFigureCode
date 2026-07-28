# =============================================================================
# Figure 3 — assembly. Loads the per-panel ggplot objects (fig3/objects/*.rds)
# and composes the final figure. Panels must be built first (fig3/panels/*.R),
# or just run fig3/run_all.R.
#
# Panel E is a reserved, undecided slot (see NOTES.md); it renders blank.
# Tagging happens HERE and only here - no panel sets labs(tag = ...).
# =============================================================================
source("fig3/00_setup.R")

P <- setNames(lapply(LETTERS[1:5],
                     function(x) readRDS(file.path(OBJ_DIR, paste0(x, ".rds")))),
              LETTERS[1:5])

coloc_counts <- P$A + P$B + plot_layout(widths = c(1, 1))
top_row      <- (coloc_counts | P$C) + plot_layout(widths = c(2.5, 1))
bottom_row   <- (P$D | P$E) + plot_layout(widths = c(1.3, 1))

Fig3 <- (top_row / bottom_row) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = SZ_TAG, face = "bold", family = FONT))

ggsave(file.path(OUTDIR, "Figure3.pdf"), Fig3,
       width = 18, height = 12, dpi = DPI, limitsize = FALSE)
# Raster preview: full physical size at lower DPI (showtext DPI matched).
showtext_opts(dpi = 100)
ggsave(file.path(OUTDIR, "Figure3.png"), Fig3,
       width = 18, height = 12, dpi = 100, limitsize = FALSE)
showtext_opts(dpi = DPI)

message("Assembled: ", file.path(OUTDIR, "Figure3.pdf"), " and Figure3.png")
