# =============================================================================
# Supplemental Figure 1 — assembly. Loads the per-panel objects
# (suppfig1/objects/*.rds) and composes the final figure.
#
# Panels A and D are each a 2-plot patchwork, so plot_annotation(tag_levels="A")
# letters SIX leaf plots (A,B on coding; C on GO; D on TF; E,F on lncRNA). That
# six-tag layout is intentional (see the design spec) and pinned by
# tests/test_suppfig1_assembly.R. Tagging happens HERE and only here.
# =============================================================================
source("suppfig1/00_setup.R")

P <- setNames(lapply(LETTERS[1:4],
                     function(x) readRDS(file.path(OBJ_DIR, paste0(x, ".rds")))),
              LETTERS[1:4])

top_row    <- P$A | (P$B + theme(plot.margin = margin(l = 2, unit = "cm")))
bottom_row <- P$C | P$D

SuppFig1 <- top_row / bottom_row +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = SZ_TAG, face = "bold", family = FONT))

ggsave(file.path(OUTDIR, "SupplementalFig1.pdf"), SuppFig1,
       width = 35, height = 23, dpi = DPI, limitsize = FALSE)
showtext_opts(dpi = 100)
ggsave(file.path(OUTDIR, "SupplementalFig1.png"), SuppFig1,
       width = 35, height = 23, dpi = 100, limitsize = FALSE)
showtext_opts(dpi = DPI)

message("Assembled: ", file.path(OUTDIR, "SupplementalFig1.pdf"), " and .png")
