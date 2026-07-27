# Extra — hemoglobin (HB*) protein-coding genes, AoU-vs-GTEx scatter. Not in the
# main figure. Inputs: aou_tpms, gtex_tpms, gencode_gtf
source("suppfig1/00_setup.R")

hb_plot <- merged_tpms() %>%
  filter(str_detect(Description, "^HB") & gene_type == "protein_coding") %>%
  ggplot(aes(x = log2(mean_TPM_AoU), y = log2(mean_TPM_GTEx))) +
  geom_point(size = 6, shape = 21, fill = "black", alpha = 0.6) +
  geom_abline(slope = 1, color = "red") +
  ylim(-5, 20) + xlim(-5, 20) +
  xlab("log2(Mean AoU TPM)") + ylab("log2(Mean GTEx TPM)") +
  ggrepel::geom_text_repel(aes(label = Description))

ggsave(file.path(OUTDIR, "SuppFig1_extra_hemoglobin.pdf"), hb_plot,
       width = 7, height = 7, dpi = DPI)
message("wrote SuppFig1_extra_hemoglobin.pdf")
