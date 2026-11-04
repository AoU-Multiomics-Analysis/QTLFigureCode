# %% CD40 colocalization and credible-set transformations
#| Start with harmonized signal results and the three SuSiE variant tables. Join only credible sets that participate in the CD40 colocalization signals.
# @include preparation/Fig4/methods.R
# @include preparation/Fig4/current_methods.R

# %% Read the harmonized CD40 signals
#| QTL_FIG4_SIGNAL can select the current signal result. Retain the CD40 gene, eQTL/pQTL/sQTL layers, and nonmissing CLPP. PIP, RCP and CLPP are posterior-based statistics.
signal_path <- prep_input_file("QTL_FIG4_SIGNAL", "input", c(
  "high_level_coloc_outputs/harmonized/harmonized_coloc.signal.tsv.gz",
  "Fig4Data/high_level_coloc_outputs/harmonized/harmonized_coloc.signal.tsv.gz"))
signal_source <- read_any(signal_path)
gene_id <- "ENSG00000101017"
cd40_signals <- signal_source %>%
  filter(gene == gene_id, layer %in% c("eQTL", "pQTL", "sQTL"), !is.na(CLPP))
message("Selected the CD40 signal results.")

# %% Read CD40 fine-mapping variants
#| Read the six required columns. Match expression and protein gene IDs exactly. Splicing IDs can contain the gene ID within the molecular-trait name. Keep credible-set and variant rows from the original results.
qtl_files <- c(eQTL = "COMB.expression_SusieMerged.annotated.ancestry_skew.tsv.gz",
               pQTL = "COMB.proteomics_SusieMerged.annotated.ancestry_skew.tsv.gz",
               sQTL = "COMB.splicing_SusieMerged..annotated.ancestry_skew.tsv.gz")
cd40_variants <- imap_dfr(qtl_files, function(relative, assay) {
  alternatives <- unique(c(relative, sub("SusieMerged[.][.]annotated", "SusieMerged.annotated", relative)))
  path <- prep_input_any("susie", alternatives)
  variants <- read_any(path, select = c("position", "variant", "cs_id", "pip", "rsid", "gene_id"))
  variants <- if (assay == "sQTL") {
    filter(variants, str_detect(.data$gene_id, fixed(.env$gene_id)))
  } else {
    filter(variants, .data$gene_id == .env$gene_id)
  }
  variants %>% mutate(layer = assay, qtl_cs = cs_id)
})
message("Read the CD40 SuSiE credible-set variants.")

# %% Join credible sets to trait and GWAS-locus records
#| Preserve the intended many-to-many expansion across trait/locus records. Calculate symmetric locus offsets within each trait and convert hg38 positions to megabases.
cd40 <- fig4_cd40_data(cd40_signals, cd40_variants, gene_id)
trait_order <- c("systemic lupus erythematosus", "multiple sclerosis", "rheumatoid arthritis",
                 "inflammatory bowel disease", "crohn's disease", "ldl cholesterol", "monocyte")
trait_order <- trait_order[trait_order %in% cd40$trait]
cd40 <- cd40 %>% mutate(assay = factor(assay, c("eQTL", "pQTL", "sQTL")),
                        credible_set = factor(credible_set, unique(credible_set)))
stopifnot(nrow(cd40) == 84L, all(cd40$gene == "ENSG00000101017"),
          all(dplyr::between(cd40$CLPP, 0, 1)), length(levels(cd40$credible_set)) <= 8L)

# %% Draw the locus display
#| The x axis gives hg38 chromosome 20 position. Color denotes assay, shape denotes the QTL credible set, and point size denotes signal-level CLPP. Calculated offsets separate multiple GWAS loci within one trait row. The dashed line marks the CD40 transcription start site.
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
