# %% Read source variant tables and refit the MAF model
#| Read the original FastENLOC variant tables and original SuSiE annotations. Apply Bayesian FDR 5% to SCP separately by assay pair. Fold every ancestry frequency to MAF, then take the maximum across ancestries. Drop missing, zero, or non-finite MAF. Fit Coloc ~ -log10(MAF) with a binomial GLM for each assay pair. Calculate 95% Wald intervals and BH-adjusted P values. The odds ratio is per tenfold decrease in MAF. No FastENLOC job is run.
# @include preparation/Fig3/methods.R
prep_result_runtime("Fig3")
variants <- coloc_variant_level()
summary <- coloc_maf_rate_fit(variants)
stopifnot(all(summary$conf.low > 0), all(summary$conf.high >= summary$conf.low),
          all(is.finite(summary$odds_ratio)))
pair_order <- c("eQTL-sQTL", "eQTL-pQTL", "sQTL-pQTL")
pair_order <- c(intersect(pair_order, summary$pair), setdiff(summary$pair, pair_order))
summary <- summary %>% mutate(pair = factor(pair, rev(pair_order)))

# %% Draw the refitted odds ratios
#| Points show the refitted unadjusted MAF odds ratio. Bars show 95% Wald intervals. This model asks whether rarer variants are more likely to colocalize. It does not adjust for fine-mapping quality and does not estimate effect-size concordance.
limits <- c(min(c(summary$conf.low, 1)) * 0.85, max(c(summary$conf.high, 1)) * 1.25)
breaks <- c(0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 8, 16)
figure <- ggplot(summary, aes(odds_ratio, pair)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), orientation = "y", width = 0, colour = "grey40") +
  geom_point(size = 2.8, colour = "#364B63") +
  scale_x_log10(limits = limits, breaks = breaks[breaks >= limits[1] & breaks <= limits[2]]) +
  labs(x = "Colocalization odds ratio per tenfold decrease in MAF", y = NULL) +
  theme_classic(base_size = 10) + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())
print(figure)

# %% Save aggregate outputs
#| Save the figure and the refitted model summaries, including raw and adjusted P values. Their interpretation remains limited by the unadjusted model.
save_plot(figure, "supplements/S23", 170, 80)
readr::write_tsv(summary, output_file("supplements/S23_model_summaries.tsv"))

message("S23 rebuilt from original result tables.")
