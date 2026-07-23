# Panel G — proportion of credible sets with a single causal variant (COMB)
source("fig1/00_setup.R")
d <- readRDS(file.path(DATA_DIR, "G.rds"))

pG <- d$summary_data %>%
  filter(group == "COMB") %>%
  left_join(d$causal_variant_count, by = c("Assay", "group")) %>%
  mutate(proportion_single = number_causal / num_cs) %>%
  mutate(Assay = factor(Assay, levels = c("Expression", "Protein", "Splicing"))) %>%
  ggplot(aes(x = Assay, y = proportion_single, fill = Assay)) +
    geom_col(show.legend = FALSE, width = 0.75, lwd = .2) +
    paletteer::scale_fill_paletteer_d("ltc::minou") +
    xlab("Assay") +
    theme(legend.position = "none",
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank()) +
    geom_text(aes(label = number_causal), vjust = 1.2, show.legend = FALSE,
              color = "black", size = 5) +
    ylab("Percent credible sets\n with single causal variant") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1))

saveRDS(pG, file.path(OBJ_DIR, "G.rds"))
message("Panel G built.")
