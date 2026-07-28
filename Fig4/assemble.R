# =============================================================================
# Figure 4 — assembly. Loads the per-panel ggplot objects (fig4/objects/*.rds)
# and composes the final figure. Panels must be built first (fig4/panels/*.R),
# or just run fig4/run_all.R.
#
# Two outputs share one layout: Figure4 (plain Panel A) and Figure4_labeled
# (Panel A with the shorthand trait annotations, from objects/A_labeled.rds).
# The labelled variant is skipped if that object has not been built.
#
# Tagging happens HERE and only here. The layout is irregular and wraps the B/C
# row in wrap_elements, so tag_levels="A" traversal is unreliable (it mislettered
# Fig3's inset); tags are assigned explicitly per panel instead.
# =============================================================================
source("fig4/00_setup.R")

P <- setNames(lapply(LETTERS[1:4],
                     function(x) readRDS(file.path(OBJ_DIR, paste0(x, ".rds")))),
              LETTERS[1:4])

compose_fig4 <- function(panel_a) {
  top_row <- ((P$B + labs(tag = "B")) | (P$C + labs(tag = "C"))) +
    plot_layout(widths = c(0.8, 1.2))
  top_row <- wrap_elements(full = top_row)   # isolate B/C from D's axis-label spacing

  right_side <- top_row / (P$D + labs(tag = "D")) +
    plot_layout(heights = c(1, 0.75))

  ((panel_a + labs(tag = "A") +
      theme(plot.margin = margin(r = 55, unit = "pt"),
            strip.text.y.left = element_text(angle = 0, hjust = 1, size = 13))) |
     right_side) +
    plot_layout(widths = c(7, 8)) &
    theme(plot.tag = element_text(size = SZ_TAG, face = "bold", family = FONT),
          plot.tag.position = c(0, 1))
}

save_fig4 <- function(fig, base) {
  ggsave(file.path(OUTDIR, paste0(base, ".pdf")), fig, width = 15, height = 10, dpi = DPI)
  showtext_opts(dpi = 100)
  ggsave(file.path(OUTDIR, paste0(base, ".png")), fig, width = 15, height = 10, dpi = 100)
  showtext_opts(dpi = DPI)
  message("Assembled: ", file.path(OUTDIR, paste0(base, ".pdf")), " and ", base, ".png")
}

Fig4 <- compose_fig4(P$A)
save_fig4(Fig4, "Figure4")

a_labeled <- file.path(OBJ_DIR, "A_labeled.rds")
if (file.exists(a_labeled)) {
  Fig4_labeled <- compose_fig4(readRDS(a_labeled))
  save_fig4(Fig4_labeled, "Figure4_labeled")
} else {
  message("Skipped Figure4_labeled: ", a_labeled, " not built ",
          "(run fig4/panels/panelA_labeled.R).")
}
