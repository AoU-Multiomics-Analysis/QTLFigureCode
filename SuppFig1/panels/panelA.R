# Panel A — protein-coding AoU-vs-GTEx TPM correlation (hex) + detection counts
# Inputs: aou_tpms, gtex_tpms, gencode_gtf
source("suppfig1/00_setup.R")

protein_coding_genes <- merged_tpms() %>%
  filter(gene_type == "protein_coding") %>%
  filter(mean_TPM_AoU > 1 | mean_TPM_GTEx > 1)

protein_coding_comparison <- protein_coding_genes %>%
  plot_gene_comparison(gene_type = "coding", hex = TRUE)

saveRDS(protein_coding_comparison, file.path(OBJ_DIR, "A.rds"))
message("Panel A built.")
