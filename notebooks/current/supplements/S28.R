# %% Read the saved expression-reduction thresholds
#| Rebuild S28 from saved distinct-variant threshold counts. The denominator is the historical set of 574 rare high-PIP eQTL variants at 453 ClinGen genes. It differs from the all-classification set in current Figure 5A/B. No variant selection or aFC model is repeated.
base <- "Supplemental Figures/01_Supplemental_drafts/Figure5_rare_eQTL_aFC_thresholds"
curve <- readr::read_tsv(paper_file(file.path(base, "threshold_curve.tsv")), show_col_types = FALSE)
landmarks <- readr::read_tsv(paper_file(file.path(base, "threshold_counts.tsv")), show_col_types = FALSE)
dots <- landmarks %>% filter(reduction_threshold_pct %in% c(0, 10, 20, 25, 40, 50, 60, 70, 75, 80, 90, 95, 100) |
                              abs(log2_afc_cutoff + 0.5) < 1e-10)
half <- landmarks %>% filter(reduction_threshold_pct == 50)
legacy <- landmarks %>% filter(abs(log2_afc_cutoff + 0.5) < 1e-10)
stopifnot(all(curve$total_assessed == 574L), all(diff(curve$variants) <= 0),
          half$variants == 45L, legacy$variants == 105L)

# %% Draw the saved thresholds
#| Expression reduction is 100 × (1 − 2^log2_aFC). A variant passes a threshold when at least one qualifying target has a strictly greater reduction. Straight lines join selected threshold points as in the current source figure; they do not estimate the intervening empirical counts.
top_breaks <- 100 * (1 - 2^c(0, -0.5, -1, -2, -log2(10)))
figure <- ggplot(dots, aes(reduction_threshold_pct, variants)) +
  geom_vline(xintercept = c(legacy$reduction_threshold_pct, 50), colour = "#A8AFB4", linewidth = 0.3, linetype = "dashed") +
  geom_line(linewidth = 0.55, colour = "#145D69") +
  geom_point(size = 1.65, colour = "#145D69") +
  geom_point(data = bind_rows(legacy, half), size = 2.5, colour = "#D55E00") +
  annotate("text", x = 1, y = 270, label = "254", hjust = 0, size = 2.8, colour = "#145D69") +
  annotate("text", x = 98, y = 270, label = "574 variants assessed", hjust = 1, size = 2.8) +
  annotate("text", x = 32, y = 137, label = "105 variants (18.3%)", hjust = 0, size = 2.8, colour = "#D55E00") +
  annotate("text", x = 54, y = 82, label = "45 variants (7.8%)\n>50% reduction", hjust = 0, size = 2.8, colour = "#D55E00") +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                     sec.axis = dup_axis(breaks = top_breaks, labels = c("0", "-0.5", "-1", "-2", "-3.32"),
                                         name = expression(log[2](aFC)~cutoff))) +
  scale_y_continuous(limits = c(0, 285), breaks = seq(0, 250, 50), expand = expansion(mult = 0)) +
  labs(x = "Expression reduction threshold (%)", y = "Rare causal eQTL variants") + theme_classic(base_size = 9)
print(figure)

# %% Save aggregate outputs
#| Save the figure and threshold counts. There is no new significance test. The 105/574 value is at log2(aFC) < −0.5, or approximately 29.3% reduction. The >50% value is 45/574.
save_plot(figure, "supplements/S28", 130, 100)
readr::write_tsv(landmarks, output_file("supplements/S28_threshold_counts.tsv"))
