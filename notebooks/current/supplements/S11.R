# %% Read the current EAS histogram
#| Rebuild S11 from the saved integer frequency table. It includes EAS and the combined analysis. No fine-mapping model is run. Traits with no credible set are absent from the source.
base <- "Supplemental Figures/01_Supplemental_drafts/Figure1_population_credible_set_counts/v2_with_EAS"
histogram <- readr::read_tsv(paper_file(file.path(base, "histogram_counts.tsv")), show_col_types = FALSE)
summary <- readr::read_tsv(paper_file(file.path(base, "population_credible_set_count_summary.tsv")), show_col_types = FALSE)
assays <- c("eQTL", "sQTL", "pQTL")
groups <- c("COMB", "AFR", "AMR", "EUR", "SAS", "EAS", "MID")
stopifnot(nrow(histogram) == 210L, nrow(summary) == 21L,
          all(histogram$number_traits >= 0), "EAS" %in% histogram$population)
check <- histogram %>% group_by(assay, population) %>%
  summarise(traits_from_bins = sum(number_traits), sets_from_bins = sum(number_traits * number_credible_sets), .groups = "drop") %>%
  left_join(summary, by = c("assay", "population"))
stopifnot(all(check$traits_from_bins == check$number_traits),
          all(check$sets_from_bins == check$total_credible_sets))

# %% Draw the complete histogram range
#| Each bar counts molecular traits with the stated number of credible sets. The 1–10 range contains all observed counts. The y scale is fixed within each assay and can differ between assays.
fills <- c(COMB = "#CF597E", AFR = "#C8D892", AMR = "#7FC385", EUR = "#30AC86", SAS = "#EBCE8D", EAS = "#ECA676", MID = "#E47D72")
ymax <- c(eQTL = 11500, sQTL = 3400, pQTL = 900)
panels <- purrr::map(assays, function(a) {
  d <- histogram %>% filter(assay == a) %>% mutate(population = factor(population, groups))
  nlabels <- summary %>% filter(assay == a) %>%
    mutate(population = factor(population, groups), label = paste0("n = ", scales::comma(number_traits)))
  ggplot(d, aes(number_credible_sets, number_traits, fill = population)) +
    geom_col(width = 0.86, colour = "grey35", linewidth = 0.15) +
    geom_text(data = nlabels, aes(x = 10.3, y = ymax[[a]] * 0.94, label = label),
              inherit.aes = FALSE, hjust = 1, vjust = 1, size = 2.5) +
    facet_grid(population ~ assay, switch = "y", labeller = labeller(population = c(COMB = "Combined", AFR = "AFR", AMR = "AMR", EUR = "EUR", SAS = "SAS", EAS = "EAS", MID = "MID"))) +
    scale_fill_manual(values = fills) +
    scale_x_continuous(breaks = 1:10, limits = c(0.4, 10.6)) +
    scale_y_continuous(limits = c(0, ymax[[a]]), labels = scales::label_comma(), expand = expansion(mult = 0)) +
    labs(x = "Credible sets per molecular trait", y = if (a == "eQTL") "Molecular traits" else NULL) +
    theme_classic(base_size = 8) +
    theme(legend.position = "none", strip.background = element_blank(), strip.placement = "outside",
          strip.text.y.left = if (a == "eQTL") element_text(angle = 0) else element_blank(),
          panel.spacing.y = grid::unit(3.5, "mm"))
})
figure <- patchwork::wrap_plots(panels, nrow = 1)
print(figure)

# %% Save aggregate outputs
#| Save the histogram and its aggregate count table. These summaries do not contain participant records.
save_plot(figure, "supplements/S11", 184, 195)
readr::write_tsv(histogram, output_file("supplements/S11_histogram.tsv"))
