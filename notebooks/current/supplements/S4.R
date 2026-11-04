# %% Read the saved population summaries
#| Rebuild the November 2 S4 display from the saved association summaries. This notebook does not rerun QTL discovery. The sample sizes are common reference-cohort counts that the original source reused across assays.
s4 <- readr::read_tsv(paper_file("Supplemental Figures/01_Supplemental_drafts/Figure1_QTL_discovery_by_population/population_QTL_discovery_data.tsv"), show_col_types = FALSE)
stopifnot(nrow(s4) == 21L, all(dplyr::between(s4$proportion_significant, 0, 1)),
          !anyDuplicated(s4[c("assay", "population")]))
s4 <- s4 %>%
  mutate(assay = factor(assay, c("eQTL", "sQTL", "pQTL")),
         population = factor(population, c("Combined", "AFR", "AMR", "EAS", "EUR", "SAS", "MID")))

# %% Draw the saved discovery proportions
#| Each point is a population and assay summary. The x axis uses the saved reference sample count. It is not an assay-specific participant count. No confidence interval or fitted curve is added.
suppressPackageStartupMessages(library(ggrepel))
fills <- s4 %>% distinct(population, colour) %>% deframe()
shapes <- s4 %>% distinct(population, marker) %>% deframe()
figure <- ggplot(s4, aes(source_reference_sample_n, proportion_significant)) +
  geom_point(aes(fill = population, shape = population), size = 2.3, colour = "#263B45", stroke = 0.35) +
  geom_text_repel(aes(label = population), seed = 71, size = 2.6,
                  min.segment.length = 0, box.padding = 0.3,
                  point.padding = 0.2, max.overlaps = Inf) +
  facet_grid(~assay) +
  scale_fill_manual(values = fills) +
  scale_shape_manual(values = shapes) +
  scale_x_log10(limits = c(280, 11000), breaks = c(300, 1000, 3000, 10000), labels = scales::label_comma()) +
  scale_y_continuous(limits = c(0, 1), labels = scales::label_percent(), breaks = seq(0, 1, 0.25)) +
  labs(x = "Reference cohort sample size (log scale)", y = "Tested genes with significant QTLs (%)") +
  theme_classic(base_size = 9) +
  theme(strip.background = element_blank(), legend.position = "none",
        panel.spacing.x = grid::unit(5, "mm"))
print(figure)

# %% Save aggregate outputs
#| Save the figure and the 21 study summaries. No participant records are read or written. The source field source_n_tested counts tested rows, not significant genes.
save_plot(figure, "supplements/S4", 184, 85)
readr::write_tsv(s4 %>% select(assay, population, proportion_significant, source_n_tested,
                              significant_count_from_saved_summary, source_reference_sample_n),
                 output_file("supplements/S4_summary.tsv"))
