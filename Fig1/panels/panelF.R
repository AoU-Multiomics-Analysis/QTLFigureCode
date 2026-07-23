# Panel F — credible-set sizes across pQTL ancestry groups (+ UKBB reference point)
source("fig1/00_setup.R")
d <- readRDS(file.path(DATA_DIR, "F.rds"))

pF <- d$pQTL_cs_sizes %>%
  PlotCSSizes(AncestryColorPalette) +
  # Black point: median credible-set size over common variants (>=1% AF, size > 1).
  geom_point(data = d$common_var_means, aes(x = group, y = mean_cs),
             shape = 21, fill = "black", size = 7, color = "black",
             show.legend = FALSE, alpha = 0.8) +
  # Diamond: external UKBB pQTL reference average.
  geom_point(data = d$ref_means, aes(x = group, y = mean_cs),
             shape = 23, size = 7, fill = "dodgerblue", color = "black") +
  geom_text_repel(inherit.aes = FALSE, data = d$ref_means,
                  aes(x = group, y = mean_cs), label = "UKBB pQTL",
                  nudge_y = 0.5, nudge_x = 0.5, color = "black", fontface = "bold",
                  segment.color = "black", segment.size = 0.6) +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank())

saveRDS(pF, file.path(OBJ_DIR, "F.rds"))
message("Panel F built.")
