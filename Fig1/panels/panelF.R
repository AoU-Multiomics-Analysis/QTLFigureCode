# Panel F — credible-set sizes across pQTL ancestry groups (+ UKBB reference point)
# Inputs: finemap_pqtl
source("fig1/00_setup.R")

pQTL_susie <- load_input("finemap_pqtl")

ref_means <- data.frame(group = "EUR", mean_cs = 9.7)   # UKBB pQTL reference average
pQTL_cs_sizes <- pQTL_susie %>%
  distinct(cs_id, cs_size, group, Assay)
common_var_means <- pQTL_susie %>%
  filter(!ALT_FREQS < .01 & !cs_size == 1) %>%
  distinct(cs_id, cs_size, group, Assay) %>%
  group_by(group) %>%
  summarize(mean_cs = median(cs_size)) %>%
  mutate(group = case_when(group == "COMB" ~ "AoU", TRUE ~ group))

pF <- pQTL_cs_sizes %>%
  PlotCSSizes(AncestryColorPalette) +
  # Black point: median credible-set size over common variants (>=1% AF, size > 1).
  geom_point(data = common_var_means, aes(x = group, y = mean_cs),
             shape = 21, fill = "black", size = 7, color = "black",
             show.legend = FALSE, alpha = 0.8) +
  # Diamond: external UKBB pQTL reference average.
  geom_point(data = ref_means, aes(x = group, y = mean_cs),
             shape = 23, size = 7, fill = "dodgerblue", color = "black") +
  geom_text_repel(inherit.aes = FALSE, data = ref_means,
                  aes(x = group, y = mean_cs), label = "UKBB pQTL",
                  nudge_y = 0.5, nudge_x = 0.5, color = "black", fontface = "bold",
                  segment.color = "black", segment.size = 0.6) +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank())

saveRDS(pF, file.path(OBJ_DIR, "F.rds"))
message("Panel F built.")
