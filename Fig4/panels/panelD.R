# Panel D — gene-set enrichment odds-ratio forest (rare vs common TWAS signal)
# Inputs: conditional_twas, clingen, drug_interactions
source("fig4/00_setup.R")

gtd <- gene_twas_data()

# Every outcome x predictor-set combination is fit; the panel plots only the
# continuous max_rare / max_common terms. The binary-predictor fits are computed
# but not plotted (see NOTES.md).
enrichment_results <- tidyr::expand_grid(outcome = names(ENRICH_OUTCOMES),
                                         model   = names(ENRICH_PREDICTORS)) %>%
  mutate(res = map2(outcome, model, ~ fit_enrich(.x, ENRICH_PREDICTORS[[.y]], gtd))) %>%
  unnest(res) %>%
  mutate(type = ENRICH_OUTCOMES[outcome]) %>%
  relocate(type, model)

plot_df <- enrichment_results %>%
  filter(term %in% c("max_rare", "max_common")) %>%
  select(type, model, term, OR = estimate, conf.low, conf.high, p.value) %>%
  mutate(term = str_remove_all(term, "max_|ever_|TRUE"),
         term = recode(term, rare = "Rare", common = "Common"),
         term = factor(term, levels = c("Rare", "Common")),
         type = factor(type, levels = c("Mendelian Gene\n(Any)",
                                        "Mendelian Gene\n(Immune/Hematological)",
                                        "Druggable gene", "Approved Drug Target")))

pal <- c(Rare = as.character(paletteer_d("futurevisions::mars"))[1], Common = "grey60")
p <- ggplot(plot_df, aes(OR, reorder(type, OR), colour = term)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high),
                  position = position_dodge(width = .55), linewidth = .8, fatten = 2.6) +
  scale_colour_manual(values = pal, name = NULL) +
  theme_classic(base_size = 12) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        strip.background = element_blank(), strip.text = element_text(face = "bold"),
        legend.position = c(0.9, 0.2), panel.spacing = unit(1.4, "lines")) +
  xlab("Odds Ratio (per unit increase in Z-score)") + ylab("")

saveRDS(p, file.path(OBJ_DIR, "D.rds"))
message("Panel D built.")
