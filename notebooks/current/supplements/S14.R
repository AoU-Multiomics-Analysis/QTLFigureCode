# %% Read the current EAS study comparisons
#| Plot the eQTL discovery comparisons. Their nine-point descriptive curves and GTEx reference envelopes are reused. This notebook does not rerun fine-mapping, count participants, or fit a new discovery curve.
base <- "Supplemental Figures/01_Supplemental_drafts/Figure1_eQTL_discovery/v2_with_EAS"
saved <- purrr::map(c("A", "B"), function(id) readRDS(paper_file(file.path(base, paste0("source_", id), "plot.rds"))))
purrr::walk(saved, function(x) {
  stopifnot(inherits(x$plot, "ggplot"), nrow(x$data) == 9L,
            "EAS" %in% x$data$population, nrow(x$gtex) == 50L)
})

# %% Restore the two discovery panels
#| Panel A uses all available targets. Panel B uses 10,881 protein-coding genes shared across the three studies. The curves describe study-level discovery and do not provide a causal estimate of sample-size effects. GTEx tissues form the saved reference envelope and are excluded from the nine-point fit.
figure <- (saved[[1]]$plot + labs(title = NULL, subtitle = NULL)) +
  (saved[[2]]$plot + labs(title = NULL, subtitle = NULL)) +
  patchwork::plot_layout(nrow = 1) + patchwork::plot_annotation(tag_levels = "A")
print(figure)

# %% Save aggregate outputs
#| Save the display and study summaries. No participant identifiers are accessed or exported.
save_plot(figure, "supplements/S14", 184, 85)
readr::write_tsv(purrr::map2_dfr(saved, c("A", "B"), ~mutate(.x$data, panel = .y)),
                 output_file("supplements/S14_study_summaries.tsv"))
