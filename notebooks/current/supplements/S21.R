# %% Read the saved gene calls
#| Rebuild S21 from saved gene-level colocalization calls. The population is all genes tested for each assay pair. It is not restricted to genes fine-mapped in all three assays. No new colocalization analysis is run.
genes <- readRDS(source_file("Fig3/data/coloc_gene_level.rds"))
stopifnot(all(c("type", "Coloc") %in% names(genes)), !anyNA(genes$Coloc))
summary <- genes %>% group_by(assay_pair = type) %>%
  summarise(n_tested = n(), n_coloc = sum(Coloc), proportion = mean(Coloc), .groups = "drop")

# %% Draw the gene-level proportions
#| Each bar is the proportion of tested genes with the saved Bayesian FDR colocalization call. A gene can contain several credible sets. These bars therefore have a different denominator from S22.
figure <- ggplot(summary, aes(assay_pair, proportion)) +
  geom_col(fill = "#364B63", width = 0.65) +
  geom_text(aes(label = scales::percent(proportion, accuracy = 0.1)), vjust = -0.4, size = 3) +
  scale_y_continuous(labels = scales::label_percent(), expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Tested genes with colocalization") + theme_classic(base_size = 10) +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank())
print(figure)

# %% Save aggregate outputs
#| Save the plot and the numerator and denominator for each assay pair. The figure has no new hypothesis test or confidence interval.
save_plot(figure, "supplements/S21", 150, 100)
readr::write_tsv(summary, output_file("supplements/S21_gene_summary.tsv"))
