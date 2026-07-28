# Panel A — gene-level colocalization: percent of genes by the number of
# colocalized molecular layers, among genes fine-mapped in all three assays.
# Inputs: enloc_gene, susie_eqtl, susie_sqtl, susie_pqtl
source("fig3/00_setup.R")

genes_all_three_assays <- susie_finemapping()$genes_all_three_assays

layer_counts <- coloc_gene_level() %>%
  filter(Gene %in% genes_all_three_assays) %>%
  group_by(Gene) %>%
  summarize(NumColocs = sum(Coloc), .groups = "drop") %>%
  dplyr::count(NumColocs) %>%
  mutate(prop = n / sum(n))

# Computed, never hardcoded: the original pasted "86.2%" into this annotation
# while computing the same value 130 lines earlier.
shared_pct <- layer_counts %>%
  filter(NumColocs > 0) %>%
  summarize(p = sum(prop)) %>%
  pull(p)

message(sprintf("Panel A: %.1f%% of genes share a signal across >=1 layer",
                100 * shared_pct))

pA <- layer_counts %>%
  mutate(NumColocs = factor(NumColocs, levels = 0:3)) %>%
  ggplot(aes(x = NumColocs, y = prop, fill = NumColocs)) +
    geom_col(width = 0.75) +
    scale_fill_manual(values = LAYER_FILL, guide = "none") +
    annotate("segment", x = 1.5, xend = 4.5, y = 0.455, yend = 0.455,
             color = "grey50", linewidth = 0.25) +
    annotate("text", x = 3, y = 0.465,
             label = sprintf("%.1f%% shared across ≥1 molecular layer",
                             100 * shared_pct),
             size = SZ_ANNOT, color = "grey20") +
    # coord_cartesian, NOT scale_y_continuous(limits = ): the latter REMOVES
    # out-of-range bars instead of clipping them, which is how panel B lost its
    # tallest bar entirely once the underlying data shifted (see NOTES.md). The
    # rendered result is identical here -- the largest share is 38.7% -- but a
    # future data change can no longer delete a bar silently.
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       breaks = seq(0, 0.5, by = 0.1)) +
    coord_cartesian(ylim = c(0, 0.55)) +
    labs(x = "Number of Colocalized Molecular Layers", y = "Percent of Genes",
         title = "Gene-level Colocalization") +
    theme(axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(hjust = 0.5, size = SZ_TITLE, face = "bold"),
          plot.title.position = "plot")

saveRDS(pA, file.path(OBJ_DIR, "A.rds"))
message("Panel A built.")
