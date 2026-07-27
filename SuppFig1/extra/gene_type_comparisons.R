# Extra — snRNA and miRNA AoU-vs-GTEx comparisons (not in the main figure).
# Inputs: aou_tpms, gtex_tpms, gencode_gtf
source("suppfig1/00_setup.R")

merged <- merged_tpms()

for (gt in c("snRNA", "miRNA")) {
  p <- merged %>%
    filter(gene_type == gt) %>%
    filter(mean_TPM_AoU > 0.1 | mean_TPM_GTEx > 0.1) %>%
    plot_gene_comparison(gene_type = gt, hex = FALSE)
  ggsave(file.path(OUTDIR, paste0("SuppFig1_extra_", gt, ".pdf")),
         p, width = 10, height = 5, dpi = DPI)
  message("wrote SuppFig1_extra_", gt, ".pdf")
}
