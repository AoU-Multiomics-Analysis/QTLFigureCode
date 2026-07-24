# Extra — number of genes detected (TPM > 0.1) per gene type in AoU. Not in the
# main figure. Inputs: aou_tpms, gencode_gtf
source("suppfig1/00_setup.R")

small_rnas_list <- c("misc_RNA", "sRNA", "vault_RNA", "scaRNA", "rRNA", "ribozyme")

# AoU per-gene means over ALL AoU-quantified genes (matches the original
# AoU_mean_TPMs_v45). Deliberately NOT merged_tpms(), which inner-joins AoU to
# GTEx and would narrow the gene universe to the AoU-GTEx intersection.
aou_means <- as_tibble(data.table::fread(input_path("aou_tpms"), skip = 2)) %>%
  compute_mean_TPM(gene_annotations())

gene_detection_summary <- aou_means %>%
  filter(mean_TPM > 0.1) %>%
  dplyr::count(gene_type) %>%
  mutate(gene_type = case_when(gene_type == "protein_coding" ~ "Coding", TRUE ~ gene_type)) %>%
  filter(!str_detect(gene_type, "IG|TR|pseudo|artifact|Mt|TEC|rRNA"))

other_small_rnas <- gene_detection_summary %>%
  filter(gene_type %in% small_rnas_list) %>%
  summarize(gene_type = "other small RNAs", n = sum(n))

detection_plot <- gene_detection_summary %>%
  bind_rows(other_small_rnas) %>%
  filter(!gene_type %in% small_rnas_list) %>%
  ggplot(aes(x = gene_type, y = n, fill = gene_type, color = gene_type)) +
  geom_col(show.legend = FALSE, alpha = 0.7) +
  ylab("Number detected \n (TPM > 0.1)") +
  scale_y_log10() +
  scale_fill_carto_d(palette = "Safe") + scale_color_carto_d(palette = "Safe") +
  geom_text(aes(label = n), vjust = -0.25, color = "black") +
  xlab("Type of gene")

ggsave(file.path(OUTDIR, "SuppFig1_extra_gene_detection.pdf"), detection_plot,
       width = 12, height = 7, dpi = DPI)
message("wrote SuppFig1_extra_gene_detection.pdf")
