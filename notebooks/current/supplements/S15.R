# %% Define annotation expansion and logistic models
#| Keep every fine-mapped variant in the background, including variants with no VEP consequence. Expand consequences into binary flags. For each flag that varies, fit annotation membership against PIP with a binomial logistic model.
# @include notebooks/current/analysis/S15_helpers.R

# %% Fit the three source annotation tables
#| Read the original combined eQTL, sQTL and pQTL annotation tables. Keep the PIP slope from each model. Adjust P values across all annotation and assay slopes with the Benjamini–Hochberg method before selecting plot annotations. Use estimate ± 1.96 × standard error for the Wald interval.
annotation_inputs <- read_tsv(repo_file("Fig1/inputs.tsv"),
  col_types = cols(.default = col_character()), show_col_types = FALSE) %>% filter(key == "annot")
stopifnot(nrow(annotation_inputs) == 3L, setequal(annotation_inputs$Assay, c("eQTL", "sQTL", "pQTL")))
fits <- map2_dfr(annotation_inputs$local_subpath, annotation_inputs$Assay, function(path, assay) {
  message("Fit annotation enrichment: ", assay)
  alternatives <- unique(c(path, sub("splicing_SusieMerged[.]annotated", "splicing_SusieMerged..annotated", path)))
  GlmAnnotationEnrichment(prep_input_any("input", alternatives)) %>% mutate(type = assay)
})
model_summaries <- s15_model_summaries(fits)
enrichment <- s15_plot_data(model_summaries)
stopifnot(all(is.finite(enrichment$odds_ratio)), all(enrichment$lower > 0),
  all(enrichment$upper >= enrichment$lower), !any(enrichment$separated))
message("Fitted source annotation models and calculated adjusted P values and Wald intervals.")

# %% Draw the annotation odds ratios
#| Points show odds ratios for PIP = 1 versus PIP = 0. Bars show the 95% Wald intervals. The dashed line marks an odds ratio of one. Opacity uses the adjusted P value.
dodge <- position_dodge(width = 0.6)
figure <- ggplot(enrichment, aes(odds_ratio, reorder(annotation_clean, odds_ratio), colour = type, alpha = significant)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.35) +
  geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = 0.15, position = dodge) +
  geom_point(size = 2, position = dodge) +
  facet_wrap(~functional_group, scales = "free_y", nrow = 1) +
  scale_x_log10(labels = scales::label_number()) +
  scale_colour_manual(values = c(eQTL = "#0072B2", sQTL = "#009E73", pQTL = "#D55E00"), name = "Assay") +
  scale_alpha_manual(values = c(`FALSE` = 0.3, `TRUE` = 1),
    labels = c(`FALSE` = "Adjusted P ≥ 0.05", `TRUE` = "Adjusted P < 0.05"), name = NULL) +
  labs(x = "Annotation enrichment odds ratio (PIP 1 vs 0)", y = NULL) +
  theme_classic(base_size = 9) +
  theme(strip.background = element_blank(), legend.position = "bottom",
    axis.line.y = element_blank(), axis.ticks.y = element_blank())
if (n_distinct(enrichment$significant) == 1L) figure <- figure + guides(alpha = "none")
print(figure)

# %% Save aggregate outputs
#| Save the plot and the fitted model summaries. Fine-mapping and variant annotation are source-result boundaries.
save_plot(figure, "supplements/S15", 184, 120)
readr::write_tsv(enrichment %>% select(type, functional_group, annotation_clean, estimate, std.error,
  conf.low, conf.high, odds_ratio, lower, upper, p.value, padj),
  output_file("supplements/S15_model_summaries.tsv"))
readr::write_tsv(model_summaries, output_file("supplements/S15_all_annotation_models.tsv"))
