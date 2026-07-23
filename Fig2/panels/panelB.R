# Panel B — population-enrichment distributions (all fine-mapped variants)
source("fig2/00_setup.R")
d <- readRDS(file.path(DATA_DIR, "B.rds"))
pd <- position_dodge(width = 0.8)

pB <- d$skew_df %>%
  ggplot(aes(pop, odds_ratio, fill = pop, color = pop,
             group = interaction(pop, rare_lab))) +
    geom_hline(yintercept = 5, linetype = "dashed", linewidth = 0.5, color = "grey40") +
    geom_quasirandom(width = 0.15, dodge.width = 0.8, size = 1, alpha = 0.7, show.legend = FALSE) +
    geom_violin(width = 0.9, alpha = 0.3, position = pd, color = NA, trim = FALSE, show.legend = FALSE) +
    stat_summary(fun = median, geom = "crossbar", linewidth = 0.5, width = 0.1,
                 alpha = 0.6, position = pd, color = "black", show.legend = FALSE) +
    scale_color_carto_d(palette = "Temps") +
    scale_fill_carto_d(palette = "Temps") +
    geom_text(data = d$skew_counts, aes(pop, 9000, label = scales::comma(n)),
              inherit.aes = FALSE, size = ANN_LG, color = "grey40", fontface = "bold") +
    geom_segment(data = d$pair_line_df, aes(x = x, xend = xend, y = y, yend = yend),
                 inherit.aes = FALSE, color = "grey60", linewidth = 1) +
    geom_text(data = d$label_df, aes(x = x, y = y, label = rare_lab),
              inherit.aes = FALSE, size = ANN_SM, color = "black", angle = 30) +
    annotate("segment", x = 0.85, xend = 7.15, y = 18000, yend = 18000,
             linewidth = 0.5, color = "grey40") +
    annotate("text", x = 4, y = 26000, label = "Ancestry-enriched variants",
             size = ANN_LG, fontface = "bold", color = "grey40") +
    scale_y_continuous(trans = "log10", breaks = c(1, 10, 100, 1000, 10000),
                       labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    labs(x = "Primary Population", y = "Causal variant population enrichment") +
    theme(axis.line.x = element_blank(), axis.ticks.x = element_blank())

saveRDS(pB, file.path(OBJ_DIR, "B.rds"))
message("Panel B built.")
