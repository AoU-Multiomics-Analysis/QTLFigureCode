# Panel B — GO enrichment of cohort-specific genes (AoU vs GTEx)
# Inputs: aou_tpms, gtex_tpms, gencode_gtf (via go_enrichment)
source("suppfig1/00_setup.R")

AoU_GTEx_enrichment_plot <- go_enrichment() %>%
  # Original: arrange(Adjusted.P.value) %>% filter(row_number() < 10) keeps
  # rows 1-9, i.e. the top 9 terms per cohort (not 10).
  format_go_plot_data(top_n = 9) %>%
  ggplot(aes(x = neglog10_p, y = reorder(Term, neglog10_p))) +
  geom_col() +
  facet_wrap(~type, scales = "free_y") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 20)) +
  xlab("-log10(FDR)") + ylab("") +
  geom_vline(xintercept = 1.3) +
  theme(panel.spacing = unit(1, "cm"), strip.text = element_text(size = 20))

saveRDS(AoU_GTEx_enrichment_plot, file.path(OBJ_DIR, "B.rds"))
message("Panel B built.")
