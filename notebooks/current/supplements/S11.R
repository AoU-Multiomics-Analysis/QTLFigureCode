# %% Source fine-mapping inputs
#| Count credible sets from variant rows for all three assays and seven analysis groups.
# @include notebooks/current/analysis/Fig1_helpers.R
fig1_groups <- c("COMB", "AFR", "AMR", "EUR", "SAS", "EAS", "MID")
# @include notebooks/current/analysis/Fig1_finemapping.R

# %% Count credible sets per molecular trait
histogram <- traits %>% count(assay, population, number_credible_sets, name = "number_traits") %>%
  complete(assay = factor(assays, levels = assays), population = factor(groups, levels = groups),
           number_credible_sets = 1:10, fill = list(number_traits = 0L)) %>%
  group_by(assay, population) %>% mutate(total_traits = sum(number_traits),
    percent_traits = 100 * number_traits / total_traits) %>% ungroup()
summary <- traits %>% group_by(assay, population) %>%
  summarise(number_traits = n(), total_credible_sets = sum(number_credible_sets),
    mean_credible_sets = mean(number_credible_sets), median_credible_sets = median(number_credible_sets),
    minimum = min(number_credible_sets), maximum = max(number_credible_sets),
    percent_one_set = 100 * mean(number_credible_sets == 1L), .groups = "drop")
stopifnot(nrow(histogram) == 210L, nrow(summary) == 21L, max(traits$number_credible_sets) == 10L)
message("Counted all source credible sets per trait and 210 histogram bins.")

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
