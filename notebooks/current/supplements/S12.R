# %% Read the saved EAS display
#| Restore S12 from the current saved plot data. No fine-mapping fit or standard-error calculation is repeated. The source contains the size distributions, the UKBB pQTL reference, and the saved AFR contrasts.
saved <- readRDS(paper_file("Supplemental Figures/01_Supplemental_drafts/Figure1_population_credible_set_sizes/v5_with_EAS/figure_display.rds"))
stopifnot(all(c("A", "B", "differences") %in% names(saved)),
          inherits(saved$A, "ggplot"), inherits(saved$B, "ggplot"),
          "EAS" %in% as.character(saved$A$data$population))

# %% Restore the two panels
#| Panel A retains the saved boxplots on the log10 size axis. Panel B shows mean size minus the AFR mean. Its bars are one gene-cluster standard error and include covariance with AFR. These are descriptive comparisons between detected signals.
panel_a <- saved$A + labs(title = NULL, subtitle = NULL)
panel_b <- saved$B + labs(title = NULL, subtitle = NULL)
figure <- panel_a / panel_b + patchwork::plot_layout(heights = c(124, 63)) +
  patchwork::plot_annotation(tag_levels = "A")
print(figure)

# %% Save aggregate outputs
#| Save the plot, exact size frequencies, and saved mean contrasts. The frequency table removes credible-set identifiers. The UKBB reference marker remains in panel A.
save_plot(figure, "supplements/S12", 184, 195)
readr::write_tsv(saved$A$data %>% count(assay, population, cs_size, name = "credible_sets"),
                 output_file("supplements/S12_size_frequencies.tsv"))
readr::write_tsv(saved$differences, output_file("supplements/S12_AFR_contrasts.tsv"))
