# %% Read the saved credible-set pairing counts
#| Rebuild S22 from saved pair-level counts. A record is a tested credible-set pairing, not a unique gene or a single-assay credible set. The source calls use Bayesian FDR 5%. No new colocalization analysis is run.
rates <- readRDS(source_file("Fig3/data/coloc_cs_rate.rds"))
stopifnot(all(rates$n_tested >= rates$n_coloc), all(rates$n_tested > 0),
          all(abs(rates$rate - rates$n_coloc / rates$n_tested) < 1e-12))

# %% Calculate the pooled rate
#| Pool counts across assay pairs. Do not average their rates. Retain the saved pair intervals and calculate the pooled 95% Wilson interval with the formula used by the current original plot builder. The saved pair intervals also match that formula.
overall <- rates %>% summarise(pair = "Overall", n_tested = sum(n_tested), n_coloc = sum(n_coloc)) %>%
  mutate(rate = n_coloc / n_tested)
z <- qnorm(0.975)
wilson <- function(n, p) {
  denominator <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denominator
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denominator
  tibble(conf.low = pmax(0, centre - half), conf.high = pmin(1, centre + half))
}
pair_check <- wilson(rates$n_tested, rates$rate)
stopifnot(max(abs(pair_check$conf.low - rates$conf.low)) < 1e-12,
          max(abs(pair_check$conf.high - rates$conf.high)) < 1e-12)
overall <- bind_cols(overall, wilson(overall$n_tested, overall$rate))
summary <- bind_rows(arrange(rates, desc(rate)), overall) %>%
  mutate(pair = factor(pair, levels = pair),
         label = paste0(scales::percent(rate, accuracy = 0.1), "\n", scales::comma(n_coloc), " / ", scales::comma(n_tested)))

# %% Draw and save the rates
#| Error bars show 95% Wilson intervals. The gray Overall bar uses the summed numerator and denominator. These intervals describe saved binary pairing calls and do not include uncertainty from fine-mapping or the colocalization model.
fills <- setNames(c("#0072B2", "#009E73", "#D55E00", "#A6A6A6"), levels(summary$pair))
figure <- ggplot(summary, aes(pair, rate, fill = pair)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.12, linewidth = 0.4) +
  geom_text(aes(y = conf.high, label = label), vjust = -0.4, size = 2.8, lineheight = 0.95) +
  scale_fill_manual(values = fills, guide = "none") +
  scale_y_continuous(labels = scales::label_percent(), expand = expansion(mult = c(0, 0.28))) +
  labs(x = NULL, y = "Credible-set pairings with colocalization") + theme_classic(base_size = 9) +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank())
print(figure)
save_plot(figure, "supplements/S22", 160, 110)
readr::write_tsv(summary %>% select(-label), output_file("supplements/S22_pairing_summary.tsv"))
