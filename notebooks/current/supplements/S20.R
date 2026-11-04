# %% Read the saved posterior results
#| Rebuild S20 from the saved FastENLOC gene, credible-set, and variant results. This notebook recalculates density displays and threshold counts from those results. It does not rerun colocalization or fine-mapping. GLCP, RCP, and SCP are posterior probabilities, not P values.
specs <- tibble(level = c("GLCP", "RCP", "SCP"),
                file = c("coloc_gene_level", "coloc_cs_level", "coloc_variant_level"))
posterior <- purrr::map2_dfr(specs$level, specs$file, function(level, file) {
  d <- readRDS(source_file(paste0("Fig3/data/", file, ".rds")))
  d %>% transmute(assay_pair = type, posterior = .data[[level]], level = level)
}) %>% mutate(level = factor(level, c("GLCP", "RCP", "SCP")))
stopifnot(all(is.finite(posterior$posterior)), all(posterior$posterior >= 0))
source_bounds <- posterior %>% group_by(level, assay_pair) %>%
  summarise(records = n(), above_one = sum(posterior > 1), maximum = max(posterior), .groups = "drop")

# %% Calculate the display summaries
#| Count records with posterior strictly greater than 0.01, 0.1, 0.5, and 0.9. The unit differs between levels: genes, credible-set pairings, or variants. The saved GLCP/RCP data include 17 values slightly above one, with maxima 1.008 and 1.002. Preserve these source values, as the original supplemental code does. Save a bounds audit rather than silently clipping them.
thresholds <- purrr::map_dfr(c(0.01, 0.1, 0.5, 0.9), function(cutoff) {
  posterior %>% group_by(level, assay_pair) %>%
    summarise(n_tested = n(), n_coloc = sum(posterior > cutoff, na.rm = TRUE), .groups = "drop") %>%
    mutate(threshold = cutoff)
})
density_plot <- ggplot(posterior, aes(posterior, fill = assay_pair)) +
  geom_density(alpha = 0.7, na.rm = TRUE) + facet_wrap(~level, nrow = 1, scales = "free_y") +
  labs(x = "Posterior probability", y = "Density", fill = "Assay pair") +
  theme_classic(base_size = 9) + theme(strip.background = element_blank(), legend.position = "bottom")
discovery_plot <- ggplot(thresholds, aes(threshold, n_coloc, colour = assay_pair)) +
  geom_point(size = 1.5) + geom_line() + facet_wrap(~level, nrow = 1) +
  scale_y_log10(labels = scales::label_comma()) +
  labs(x = "Posterior threshold", y = "Colocalization records (log scale)", colour = "Assay pair") +
  theme_classic(base_size = 9) + theme(strip.background = element_blank(), legend.position = "none")

# %% Display and save aggregate results
#| The two rows match the scope of the original supplemental QC code. Save only the figure and threshold summaries. No variant or participant rows are exported.
figure <- density_plot / discovery_plot + patchwork::plot_annotation(tag_levels = "A")
print(figure)
save_plot(figure, "supplements/S20", 184, 160)
readr::write_tsv(thresholds, output_file("supplements/S20_threshold_counts.tsv"))
readr::write_tsv(source_bounds, output_file("supplements/S20_source_bounds.tsv"))
