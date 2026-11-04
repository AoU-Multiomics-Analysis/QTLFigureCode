# %% Read the saved trait ranks
#| Rebuild S30 from four saved signed-effect percentile ranks. No association analysis, clinical classification, or ranking is repeated. The input contains four gene–trait summaries and no participant records.
panels <- readr::read_tsv(paper_file("Supplemental Figures/01_Supplemental_drafts/Figure5_trait_prioritization/trait_prioritization_data.tsv"), show_col_types = FALSE) %>%
  mutate(gene = factor(gene, c("F10", "ANGPTL3", "SERPINF2", "IFT80")))
stopifnot(nrow(panels) == 4L, !anyNA(panels), all(dplyr::between(panels$percentile, 0, 100)),
          all(abs(panels$x - if_else(panels$tail == "Lower tail", panels$percentile / 5, (panels$percentile - 95) / 5)) < 1e-12))
ticks <- tidyr::crossing(gene = panels$gene, x = seq(0, 1, 0.2)) %>%
  mutate(label = if_else(gene == "ANGPTL3", x * 5, 95 + x * 5))

# %% Draw the four selected examples
#| F10, SERPINF2, and IFT80 show the 95th–100th percentile range. ANGPTL3 shows the 0th–5th range. The upper-tail percentage is 100 minus the percentile; the lower-tail percentage is the percentile. These directional ranks are not association P values or individual disease-risk estimates.
figure <- ggplot() +
  geom_text(data = panels, aes(0.5, 1.24, label = gene), size = 3.5, fontface = "bold") +
  geom_text(data = panels, aes(0.5, 1.04, label = trait), size = 2.9) +
  geom_text(data = panels, aes(0.5, 0.85, label = ancestry), size = 2.5, colour = "#52676D") +
  annotate("segment", x = 0, xend = 1, y = 0.72, yend = 0.72, colour = "#C7D0D3", linewidth = 0.3) +
  geom_text(data = panels, aes(0.5, 0.49, label = label), size = 3.2, fontface = "bold", colour = "#263B45") +
  geom_text(data = panels, aes(0.5, 0.31, label = rank_label), size = 2.5, colour = "#52676D") +
  annotate("segment", x = 0, xend = 1, y = 0.06, yend = 0.06, colour = "#C6C9CC", linewidth = 0.7) +
  geom_segment(data = panels, aes(x = x, xend = tail_end, colour = class), y = 0.06, yend = 0.06, linewidth = 0.85, show.legend = FALSE) +
  geom_point(data = panels, aes(x, 0.06, colour = class, shape = class), size = 2.7) +
  geom_segment(data = ticks, aes(x = x, xend = x), y = -0.09, yend = -0.13, colour = "#73777B", linewidth = 0.25) +
  geom_text(data = ticks, aes(x, -0.24, label = label), size = 2.5, colour = "#52676D") +
  geom_text(data = panels, aes(0.5, -0.44, label = tail), size = 2.5, colour = "#52676D") +
  facet_wrap(vars(gene), ncol = 2) +
  scale_colour_manual(values = c(VUS = "#E88471", Unclassified = "#009392"), name = NULL) +
  scale_shape_manual(values = c(VUS = 17, Unclassified = 15), name = NULL) +
  labs(x = "Genome-wide signed-effect percentile") +
  coord_cartesian(xlim = c(-0.07, 1.07), ylim = c(-0.51, 1.36), expand = FALSE, clip = "off") +
  theme_void(base_size = 9) +
  theme(strip.text = element_blank(), panel.spacing.x = grid::unit(7, "mm"), panel.spacing.y = grid::unit(5, "mm"),
        axis.title.x = element_text(size = 9, margin = margin(t = 4)), legend.position = "bottom")
print(figure)

# %% Save aggregate outputs
#| Save the figure and the four source summaries. VUS and unclassified labels retain the source classification and are not a new clinical assessment.
save_plot(figure, "supplements/S30", 184, 125)
readr::write_tsv(panels, output_file("supplements/S30_trait_ranks.tsv"))
