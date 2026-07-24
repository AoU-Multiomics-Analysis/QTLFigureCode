# Extra — rRNA mean-TPM histogram (AoU vs GTEx). Not in the main figure.
# Inputs: aou_tpms, gtex_tpms, gencode_gtf
source("suppfig1/00_setup.R")

rRNA_hist <- merged_tpms() %>%
  filter(gene_type == "rRNA") %>%
  select(Description, mean_TPM_AoU, mean_TPM_GTEx) %>%
  pivot_longer(!Description) %>%
  mutate(name = str_remove(name, "mean_TPM_")) %>%
  ggplot(aes(x = value, fill = name)) +
  geom_histogram(position = "identity", alpha = 0.7) +
  xlab("Mean rRNA TPM") + ylab("Number genes") +
  scale_fill_carto_d(palette = "Safe")

ggsave(file.path(OUTDIR, "SuppFig1_extra_rRNA.pdf"), rRNA_hist,
       width = 7, height = 5, dpi = DPI)
message("wrote SuppFig1_extra_rRNA.pdf")
