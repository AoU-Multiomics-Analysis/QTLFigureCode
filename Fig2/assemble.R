# =============================================================================
# Figure 2 — assembly. Loads the per-panel ggplot objects (fig2/objects/*.rds)
# and composes the final figure (Draft2 layout). Panels must be built first
# (fig2/01_prep_data.R then fig2/panels/panel[A-G].R), or just run fig2/run_all.R.
# =============================================================================
source("fig2/00_setup.R")

P <- setNames(lapply(LETTERS[1:7], function(x) readRDS(file.path(OBJ_DIR, paste0(x, ".rds")))),
              LETTERS[1:7])

# E over F form the QTL-burden column.
QTLBurden <- P$E / P$F + plot_layout(heights = c(1, 3))

top_row <- P$A + (P$B + theme(plot.margin = margin(l = 100))) + P$C +
  plot_layout(widths = c(1, 2.5, 1))

Fig2 <- top_row / (P$D | QTLBurden | P$G) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = SZ_TAG, face = "bold", family = FONT))

ggsave(file.path(OUTDIR, "Figure2.pdf"), Fig2, width = 40, height = 24, dpi = DPI, limitsize = FALSE)
# Raster preview: full physical size at lower DPI (showtext DPI matched).
showtext_opts(dpi = 100)
ggsave(file.path(OUTDIR, "Figure2.png"), Fig2, width = 40, height = 24, dpi = 100, limitsize = FALSE)
showtext_opts(dpi = DPI)

message("Assembled: ", file.path(OUTDIR, "Figure2.pdf"), " and Figure2.png")
