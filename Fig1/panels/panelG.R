# Panel G — proportion of credible sets with a single causal variant (COMB)
# Inputs: finemap_eqtl, finemap_pqtl, finemap_sqtl
source("fig1/00_setup.R")

finemapping_cs <- bind_rows(
  load_input("finemap_eqtl") %>% select(group, pip, Assay, cs_id),
  load_input("finemap_pqtl") %>% select(group, pip, Assay, cs_id),
  load_input("finemap_sqtl") %>% select(group, pip, Assay, cs_id)
)

summary_data <- finemapping_cs %>%
  distinct(cs_id, group, Assay) %>%
  dplyr::count(Assay, group) %>%
  dplyr::rename("num_cs" = "n")
causal_variant_count <- finemapping_cs %>%
  filter(pip > 0.9) %>%
  distinct(cs_id, group, Assay) %>%
  dplyr::count(Assay, group) %>%
  dplyr::rename("number_causal" = "n")

pG <- summary_data %>%
  filter(group == "COMB") %>%
  left_join(causal_variant_count, by = c("Assay", "group")) %>%
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
