# =============================================================================
# Figure 4 — assembly. Loads the per-panel ggplot objects (fig4/objects/*.rds)
# and composes the final figure. Panels must be built first (fig4/panels/*.R),
# or just run fig4/run_all.R.
#
# Tagging happens HERE and only here. The layout is irregular and wraps the B/C
# row in wrap_elements, so tag_levels="A" traversal is unreliable (it mislettered
# Fig3's inset); tags are assigned explicitly per panel instead.
# =============================================================================
source("fig4/00_setup.R")

P <- setNames(lapply(LETTERS[1:4],
                     function(x) readRDS(file.path(OBJ_DIR, paste0(x, ".rds")))),
              LETTERS[1:4])

top_row <- ((P$B + labs(tag = "B")) | (P$C + labs(tag = "C"))) +
  plot_layout(widths = c(0.8, 1.2))
top_row <- wrap_elements(full = top_row)   # isolate B/C from D's axis-label spacing

right_side <- top_row / (P$D + labs(tag = "D")) +
  plot_layout(heights = c(1, 0.75))

Fig4 <- ((P$A + labs(tag = "A") +
            theme(plot.margin = margin(r = 55, unit = "pt"),
                  strip.text.y.left = element_text(angle = 0, hjust = 1, size = 13))) |
           right_side) +
  plot_layout(widths = c(7, 8)) &
  theme(plot.tag = element_text(size = SZ_TAG, face = "bold", family = FONT),
        plot.tag.position = c(0, 1))

ggsave(file.path(OUTDIR, "Figure4.pdf"), Fig4, width = 15, height = 10, dpi = DPI)
showtext_opts(dpi = 100)
ggsave(file.path(OUTDIR, "Figure4.png"), Fig4, width = 15, height = 10, dpi = 100)
showtext_opts(dpi = DPI)

message("Assembled: ", file.path(OUTDIR, "Figure4.pdf"), " and Figure4.png")
