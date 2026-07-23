# Panel C — proportion of genes with significant QTLs (expression | splicing + protein)
# Inputs: perm (all assays), ancestry_preds, comb_sample_list
source("fig1/00_setup.R")

sample_ids <- fread(input_path("comb_sample_list"))$research_id
ancestry_counts <- LoadAncestryCounts(input_path("ancestry_preds"), sample_ids)

# FDR is per (ancestry, assay) — i.e. within each permutation file.
QTL_data <- load_input("perm") %>%
  group_by(group, Assay) %>%
  mutate(FDR = p.adjust(pval_beta, method = "fdr")) %>%
  ungroup()

combined_count <- sum(ancestry_counts$n)
prop_sig_data <- QTL_data %>%
  SummarizeSigQTLs(ancestry_counts) %>%
  left_join(AncestryColorPalette, by = c("group" = "Group")) %>%
  mutate(n = case_when(group == "AoU" ~ combined_count, TRUE ~ n))

add_aou_label <- function(df) {
  df %>% mutate(label = case_when(group == "AoU" ~ paste0("AoU", "\n(", num_sig, ")"), TRUE ~ group))
}
expression_only     <- prop_sig_data %>% filter(Assay == "Expression") %>% add_aou_label()
splicing_proteomics <- prop_sig_data %>% filter(Assay != "Expression") %>% add_aou_label()

expression_plot <- expression_only %>%
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
    data = subset(expression_only, Assay == "Expression" & group == "AoU"),
    aes(label = label),
    size = 8, force = 1000, nudge_x = -.4, nudge_y = -.025,
    segment.color = "black", segment.size = 1.5, segment.alpha = 0.7,
    min.segment.length = 0
  ) +
  guides(shape = guide_legend(override.aes = list(size = 2.5)))

splicing_proteomics_plot <- splicing_proteomics %>%
  PlotPropQTLSig() +
  theme(strip.text = element_text(size = 14),
        axis.title.y = element_blank(),
        plot.tag = element_blank()) +
  guides(fill = "none", shape = "none")

# Combined object: expression (tagged) beside splicing+protein (tag suppressed).
pC <- expression_plot + splicing_proteomics_plot +
  plot_layout(axis_titles = "collect", widths = c(2, 1))

saveRDS(pC, file.path(OBJ_DIR, "C.rds"))
message("Panel C built.")
