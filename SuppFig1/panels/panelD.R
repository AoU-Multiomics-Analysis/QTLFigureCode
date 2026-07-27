# Panel D — lncRNA AoU-vs-GTEx TPM correlation (hex) + detection counts
# Inputs: aou_tpms, gtex_tpms, gencode_gtf
source("suppfig1/00_setup.R")

lncRNA_comparison <- merged_tpms() %>%
  filter(gene_type == "lncRNA") %>%
  filter(mean_TPM_AoU > 0.1 | mean_TPM_GTEx > 0.1) %>%
  plot_gene_comparison(gene_type = "lncRNA", hex = TRUE)

saveRDS(lncRNA_comparison, file.path(OBJ_DIR, "D.rds"))
message("Panel D built.")
