# Panel A — MAF cumulative distribution of causal variants
source("fig2/00_setup.R")
d <- readRDS(file.path(DATA_DIR, "A.rds"))

pA <- d$df %>%
  ggplot(aes(log10(gvs_max_af), color = Assay)) +
    annotate("rect", xmin = -Inf, xmax = log10(0.01), ymin = -Inf, ymax = Inf,
             alpha = 0.1, fill = "grey70") +
    stat_ecdf(linewidth = 2.5, alpha = 0.8) +
    paletteer::scale_color_paletteer_d("ltc::minou") +
    scale_x_continuous(breaks = c(-1, -2, -3), labels = c("0.1", "0.01", "0.001")) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    annotate("text", x = d$x_center, y = Inf, label = "Rare variants",
             color = "grey40", vjust = 1.5, size = ANN_LG, fontface = "bold") +
    labs(x = "MAF", y = "Cumulative % of causal\nvariants (MAF ≤ k)") +
    theme(legend.position = c(0.85, 0.2), axis.line.x = element_blank())

saveRDS(pA, file.path(OBJ_DIR, "A.rds"))
message("Panel A built.")
