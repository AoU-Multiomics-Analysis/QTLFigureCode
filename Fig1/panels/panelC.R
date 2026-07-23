# Panel C — proportion of genes with significant QTLs (expression | splicing + protein)
source("fig1/00_setup.R")
d <- readRDS(file.path(DATA_DIR, "C.rds"))

expression_plot <- d$expression_only %>%
  PlotPropQTLSig() +
  theme(
    legend.text  = element_text(size = 15),
    legend.title = element_text(size = 15),
    panel.spacing = unit(2, "lines"),
    legend.position = c(0.88, 0.22),
    legend.key.size = unit(0.2, "cm"),
    plot.margin = margin(r = 20, unit = "pt")
  ) +
  ggrepel::geom_text_repel(
    data = subset(d$expression_only, Assay == "Expression" & group == "AoU"),
    aes(label = label),
    size = 8, force = 1000, nudge_x = -.4, nudge_y = -.025,
    segment.color = "black", segment.size = 1.5, segment.alpha = 0.7,
    min.segment.length = 0
  ) +
  guides(shape = guide_legend(override.aes = list(size = 2.5)))

splicing_proteomics_plot <- d$splicing_proteomics %>%
  PlotPropQTLSig() +
  theme(strip.text = element_text(size = 14),
        axis.title.y = element_blank(),
        plot.tag = element_blank()) +
  guides(fill = "none", shape = "none")

# Combined object: expression (tagged C) beside splicing+protein (tag suppressed).
pC <- expression_plot + splicing_proteomics_plot +
  plot_layout(axis_titles = "collect", widths = c(2, 1))

saveRDS(pC, file.path(OBJ_DIR, "C.rds"))
message("Panel C built.")
