# %% Read the saved enrichment estimates
#| Rebuild S15 from the logistic-regression summaries stored in the saved Figure 1 enrichment plot. No regression is fitted here. Odds ratios compare PIP = 1 with PIP = 0 within each assay.
saved <- readRDS(source_file("Fig1/objects/H.rds"))
enrichment <- saved$data %>%
  mutate(odds_ratio = exp(estimate), lower = exp(conf.low), upper = exp(conf.high),
         type = factor(type, c("eQTL", "sQTL", "pQTL")), significant = padj < 0.05)
stopifnot(all(is.finite(enrichment$odds_ratio)), all(enrichment$lower > 0),
          all(enrichment$upper >= enrichment$lower), !any(enrichment$separated))

# %% Draw the saved odds ratios
#| Points show saved odds ratios. Bars show saved 95% Wald intervals. The log10 axis displays odds ratios rather than log-odds coefficients. The alpha key uses the saved adjusted P values.
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
  theme(strip.background = element_blank(), legend.position = "bottom", axis.line.y = element_blank(), axis.ticks.y = element_blank())
if (n_distinct(enrichment$significant) == 1L) figure <- figure + guides(alpha = "none")
print(figure)

# %% Save aggregate outputs
#| Export the figure and saved model summaries. Fine-mapping, variant annotation, and logistic-model fitting remain upstream steps.
save_plot(figure, "supplements/S15", 184, 120)
readr::write_tsv(enrichment %>% select(type, functional_group, annotation_clean, estimate,
                                      conf.low, conf.high, odds_ratio, lower, upper, p.value, padj),
                 output_file("supplements/S15_model_summaries.tsv"))
