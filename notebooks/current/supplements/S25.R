# %% Read the saved CD40 locus data
#| Rebuild S25 from the saved CD40 display data. The source already links QTL credible sets to GWAS trait and locus records. No harmonization, fine-mapping, or colocalization is repeated. PIP, RCP, and CLPP are posterior-based statistics, not association P values.
cd40 <- readr::read_tsv(source_file("Fig4/cd40_coloc_locus_plot_data.tsv"), show_col_types = FALSE)
trait_order <- c("systemic lupus erythematosus", "multiple sclerosis", "rheumatoid arthritis",
                 "inflammatory bowel disease", "crohn's disease", "ldl cholesterol", "monocyte")
trait_order <- trait_order[trait_order %in% cd40$trait]
cd40 <- cd40 %>% mutate(assay = factor(assay, c("eQTL", "pQTL", "sQTL")),
                        credible_set = factor(credible_set, unique(credible_set)))
stopifnot(nrow(cd40) > 0, all(cd40$gene == "ENSG00000101017"),
          all(dplyr::between(cd40$CLPP, 0, 1)), length(levels(cd40$credible_set)) <= 8L)

# %% Draw the locus display
#| The x axis gives hg38 chromosome 20 position. Color denotes assay, shape denotes the QTL credible set, and point size denotes signal-level CLPP. Saved offsets separate multiple GWAS loci within one trait row. The dashed line marks the CD40 transcription start site.
shapes <- setNames(c(16, 15, 18, 17, 8, 16, 15, 18)[seq_along(levels(cd40$credible_set))], levels(cd40$credible_set))
x_limits <- range(cd40$position_mb, finite = TRUE)
x_limits <- x_limits + c(-1, 1) * max(0.003, diff(x_limits) * 0.06)
figure <- ggplot(cd40, aes(position_mb, trait_y)) +
  geom_vline(xintercept = 46.118278, linetype = "dashed", linewidth = 0.35, colour = "grey50") +
  geom_point(aes(shape = credible_set, colour = assay, size = CLPP), alpha = 0.9) +
  scale_shape_manual(values = shapes, name = "QTL credible set") +
  scale_colour_manual(values = c(eQTL = "#0072B2", pQTL = "#D55E00", sQTL = "#009E73"), name = "Assay") +
  scale_size_continuous(range = c(2.5, 8.5), limits = c(0, 1), breaks = c(0.05, 0.1, 0.2, 0.4), name = "CLPP") +
  scale_x_continuous(limits = x_limits, breaks = pretty(x_limits, 6), labels = scales::label_number(accuracy = 0.001)) +
  scale_y_continuous(breaks = seq_along(trait_order), labels = rev(trait_order), expand = expansion(mult = c(0.08, 0.08))) +
  labs(x = "Chromosome 20 position (Mb)", y = NULL) +
  guides(shape = guide_legend(override.aes = list(size = 3, colour = "grey40"), nrow = 2),
         colour = guide_legend(override.aes = list(shape = 16, size = 3)), size = guide_legend(order = 4)) +
  theme_classic(base_size = 9) +
  theme(axis.text.y = element_text(size = 8), legend.text = element_text(size = 7.5),
        panel.grid.major.y = element_line(colour = "#EEEEEE", linewidth = 0.3))
print(figure)

# %% Save the figure and aggregate audit
#| Save a trait–assay summary rather than variant-level source records. Variant rows can repeat across GWAS loci and must not be interpreted as independent discoveries.
save_plot(figure, "supplements/S25", 230, 165)
readr::write_tsv(cd40 %>% group_by(trait, assay) %>%
                   summarise(variants = n_distinct(variant), credible_sets = n_distinct(qtl_cs),
                             loci = n_distinct(locus_key), maximum_CLPP = max(CLPP), .groups = "drop"),
                 output_file("supplements/S25_trait_assay_summary.tsv"))
