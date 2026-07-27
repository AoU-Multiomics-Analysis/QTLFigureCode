# Panel C — TF expression boxplot (AoU vs GTEx), log2(TPM+1) once
# Inputs: aou_tpms, gtex_tpms, gencode_gtf, tf_list
source("suppfig1/00_setup.R")

TF_genes <- tf_expression_table(merged_tpms(), tf_list())

# Reported statistics (log exactly once, matching the axis label).
tf_lm <- TF_genes %>%
  do(broom::tidy(lm(log2(mean_TPM_GTEx + 1) ~ log2(mean_TPM_AoU + 1), data = .)))
tf_ttest <- t.test(log2(TF_genes$mean_TPM_AoU + 1),
                   log2(TF_genes$mean_TPM_GTEx + 1), paired = TRUE)
message("Panel C TF lm:"); print(tf_lm)
message("Panel C TF paired t-test p = ", signif(tf_ttest$p.value, 3))

TF_gene_comparison <- TF_genes %>%
  select(Description, contains("TPM")) %>%
  pivot_longer(!Description) %>%
  mutate(name = str_remove(name, "mean_TPM_")) %>%
  ggplot(aes(x = name, y = log2(value + 1), fill = name)) +
  geom_boxplot(show.legend = FALSE) +
  ylab("TF Expression \n log2(TPM + 1)") +
  scale_fill_carto_d(palette = "Safe")

saveRDS(TF_gene_comparison, file.path(OBJ_DIR, "C.rds"))
message("Panel C built.")
