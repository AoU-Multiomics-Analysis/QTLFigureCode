# %% Read the saved MAF model estimates
#| Rebuild S23 from saved logistic-model summaries. The odds ratio is per tenfold decrease in MAF. No logistic model or colocalization analysis is fitted here.
summary <- readRDS(source_file("Fig3/data/coloc_maf_rate.rds"))
stopifnot(all(summary$conf.low > 0), all(summary$conf.high >= summary$conf.low),
          all(is.finite(summary$odds_ratio)))
pair_order <- c("eQTL-sQTL", "eQTL-pQTL", "sQTL-pQTL")
pair_order <- c(intersect(pair_order, summary$pair), setdiff(summary$pair, pair_order))
summary <- summary %>% mutate(pair = factor(pair, rev(pair_order)))

# %% Draw the saved odds ratios
#| Points show the saved unadjusted MAF odds ratio. Bars show 95% Wald intervals. This model asks whether rarer variants are more likely to colocalize. It does not adjust for fine-mapping quality and does not estimate effect-size concordance.
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
#| Save the figure and the saved model summaries, including raw and adjusted P values. Their interpretation remains limited by the unadjusted model.
save_plot(figure, "supplements/S23", 170, 80)
readr::write_tsv(summary, output_file("supplements/S23_model_summaries.tsv"))
