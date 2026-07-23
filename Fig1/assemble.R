# =============================================================================
# Figure 1 — assembly. Loads the per-panel ggplot objects (fig1/objects/*.rds)
# and composes the final figure (4-row layout from QTL_paper_Fig1.ipynb). Panels
# must be built first (fig1/01_prep_data.R then fig1/panels/panel[A-H].R), or
# just run fig1/run_all.R.
#
# Panel A is a blank slot for the study-design schematic added externally.
# =============================================================================
source("fig1/00_setup.R")

P <- setNames(lapply(LETTERS[1:8], function(x) readRDS(file.path(OBJ_DIR, paste0(x, ".rds")))),
              LETTERS[1:8])

top_row    <- P$A + (P$B + theme(plot.margin = margin(l = 60, unit = "pt")))
second_row <- (P$C + theme(plot.margin = margin(r = 60, b = 40, unit = "pt"))) |
              (P$D + theme(plot.margin = margin(b = 40, unit = "pt")))
third_row  <- P$E | (P$F + theme(plot.margin = margin(l = 20, b = 20, unit = "pt")))
bottom_row <- (P$G + theme(plot.margin = margin(r = 60, unit = "pt"))) + P$H +
              plot_layout(widths = c(1, 2.5))

Fig1 <- (top_row / second_row / third_row / bottom_row) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = SZ_TAG))

ggsave(file.path(OUTDIR, "Figure1.pdf"), Fig1, width = 25, height = 27, dpi = DPI, limitsize = FALSE)
# Raster preview: full physical size at lower DPI (showtext DPI matched).
showtext_opts(dpi = 100)
ggsave(file.path(OUTDIR, "Figure1.png"), Fig1, width = 25, height = 27, dpi = 100, limitsize = FALSE)
showtext_opts(dpi = DPI)

message("Assembled: ", file.path(OUTDIR, "Figure1.pdf"), " and Figure1.png")
