# %% Read the saved REM display
#| Restore S16 from the current saved ancestry-enrichment plot. The source already displays OTH as REM. This notebook does not rerun ancestry enrichment tests or estimate the rolling summaries again.
saved <- readRDS(paper_file("Supplemental Figures/01_Supplemental_drafts/Figure2_ancestry_enrichment_REM/Supplemental_Figure_S16_ancestry_enrichment_REM.rds"))
stopifnot(inherits(saved, "ggplot"),
          all(c("af_bin", "n_enriched", "total_variants", "prop", "x") %in% names(saved$data)),
          all(dplyr::between(saved$data$prop, 0, 1)))

# %% Restore the ancestry-frequency curves
#| Colored curves show ancestry-enriched variant proportions in saved frequency bins. The additional curve shows variants with nonzero frequency in all seven groups. Its saved rolling mean is displayed with the original unconstrained LOESS layer, which can exceed 100%. The smooth curve is a visual guide. No enrichment test is repeated.
figure <- saved + labs(title = NULL, subtitle = NULL)
print(figure)

# %% Save aggregate outputs
#| Save the figure and frequency-bin counts. These are aggregate variant counts. REM is a display label for the source OTH category, not a new population estimate.
save_plot(figure, "supplements/S16", 184, 125)
readr::write_tsv(saved$data %>% select(af_bin, gvs_max_subpop, n_enriched, total_variants, prop, x) %>%
                   mutate(gvs_max_subpop = recode(as.character(gvs_max_subpop), OTH = "REM")),
                 output_file("supplements/S16_frequency_bins.tsv"))
