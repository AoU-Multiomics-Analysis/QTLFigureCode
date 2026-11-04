# %% Read the saved EAS PIP panels
#| Restore S13 from the current density curves and summaries. This notebook does not estimate new densities or run fine-mapping. The density grid contains saved display coordinates, not variant records.
saved <- readRDS(paper_file("Supplemental Figures/01_Supplemental_drafts/Figure1_population_lead_PIP/v2_with_EAS/figure_display.rds"))
stopifnot(all(c("A", "B", "C", "summary") %in% names(saved)),
          "EAS" %in% as.character(saved$summary$population),
          all(dplyr::between(saved$A$data$pip, 0, 1)))

# %% Restore density and summary panels
#| Panel A shows the saved lead-PIP densities. Panels B and C show the mean lead PIP and the percentage with lead PIP strictly greater than 0.9. Population analyses can contain different detected credible sets.
figure <- (saved$A + labs(title = NULL, subtitle = NULL)) /
  (saved$B + labs(title = NULL, subtitle = NULL)) /
  (saved$C + labs(title = NULL, subtitle = NULL)) +
  patchwork::plot_layout(heights = c(1.3, 1, 1)) + patchwork::plot_annotation(tag_levels = "A")
print(figure)

# %% Save aggregate outputs
#| Save the figure and population summaries. Percentages use a 0–100 scale.
save_plot(figure, "supplements/S13", 184, 210)
readr::write_tsv(saved$summary, output_file("supplements/S13_summary.tsv"))
