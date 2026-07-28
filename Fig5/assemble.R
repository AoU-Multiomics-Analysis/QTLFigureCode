# =============================================================================
# Figure 5 — assembly. Loads the per-panel objects (fig5/objects/*.rds) and
# composes the final figure. Build panels first (fig5/panels/*.R), or just run
# fig5/run_all.R.
#
# TAGGING LIVES HERE (only here). The panels are a heterogeneous ggplot /
# patchwork-with-inset / cowplot mix, so the composite panels (A, D, F) are
# wrapped with wrap_elements() to each count as a single tag leaf (G is already a
# wrap_elements; B/C/E are bare ggplots). Tags are then assigned by an EXPLICIT
# sequence that matches patchwork's depth-first traversal of the row-major layout
# — A,B,C across the top, then the fine-map row (D,F) above the mechanism row
# (E,G) — so the visual panels read A-G. Verified by tests/test_fig5_assembly.R.
#
# Layout: top row = A | (B|C); bottom half = fine-map row (D|F) over mechanism
# row (E|G). D and F share a row so they render at the same size.
# =============================================================================
source("fig5/00_setup.R")

P <- setNames(lapply(LETTERS[1:7],
                     function(x) readRDS(file.path(OBJ_DIR, paste0(x, ".rds")))),
              LETTERS[1:7])

# Collapse the composite panels into single tag leaves (G is already wrapped).
A_panel <- wrap_elements(full = P$A)
D_panel <- wrap_elements(full = P$D)
F_panel <- wrap_elements(full = P$F)

combined_BC <- (P$B | P$C) + plot_layout(widths = c(0.85, 1.05))
top_row     <- A_panel | combined_BC + plot_layout(widths = c(1.5, 1))

finemap_row   <- D_panel | F_panel
mechanism_row <- P$E | P$G
# Row heights are written as inches of the 16-inch canvas (patchwork normalises
# them, so only the ratios matter, but the inches are what we reason about).
# D/F were too squat for the regional plot + inset + gene track they carry, so the
# fine-map row takes ~0.8 in from the top row and the mechanism row; the spacers
# stay wide enough to read as gutters between the three bands.
bottom_half   <- finemap_row / plot_spacer() / mechanism_row +
  plot_layout(heights = c(4.00, 0.26, 6.99))

Fig5 <- top_row / plot_spacer() / bottom_half +
  plot_layout(heights = c(4.35, 0.40, 11.25)) +
  plot_annotation(tag_levels = list(c("A", "B", "C", "D", "F", "E", "G"))) &
  theme(plot.tag = element_text(size = SZ_TAG, face = "bold", family = FONT))

ggsave(file.path(OUTDIR, "Figure5.pdf"), Fig5, width = 18, height = 16, bg = "white",
       dpi = DPI, limitsize = FALSE)
showtext_opts(dpi = 100)
ggsave(file.path(OUTDIR, "Figure5.png"), Fig5, width = 18, height = 16, bg = "white",
       dpi = 100, limitsize = FALSE)
showtext_opts(dpi = DPI)

message("Assembled: ", file.path(OUTDIR, "Figure5.pdf"), " and Figure5.png")
