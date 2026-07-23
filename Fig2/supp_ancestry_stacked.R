# =============================================================================
# Supplemental figure — composition of ancestry-enriched variants vs. allele
# frequency (expression). For each MAF bin: the fraction of variants that are
# ancestry-enriched, split by enriched ancestry (stacked areas), with the
# proportion shared across all continents overlaid (black line).
#
# Built over fine-mapped expression variants (PIP > 0.9). To use all SuSiE-tested
# variants instead, remove the `filter(pip > 0.9)` line marked below.
#
# Run from project root:  Rscript fig2/supp_ancestry_stacked.R
# =============================================================================
source("fig2/00_setup.R")

message("Loading expression input ...")
eQTL_COMB_susie <- add_derived(read_qtl("COMB.expression_SusieMerged.annotated.ancestry_skew.tsv.gz"))

af_cols <- c("gvs_all_af", "gvs_max_af", "gvs_afr_af", "gvs_amr_af",
             "gvs_eas_af", "gvs_eur_af", "gvs_mid_af", "gvs_sas_af", "gvs_oth_af")

AncestrySkewAnnotated <- eQTL_COMB_susie %>%
  filter(pip > 0.9) %>%                        # fine-mapped variants only
  mutate(across(all_of(af_cols), as.numeric)) %>%
  mutate(
    MAF     = if_else(gvs_all_af > 0.5, 1 - gvs_all_af, gvs_all_af),
    log_MAF = -log10(MAF),
    padj    = p.adjust(p_value, method = "fdr"),
    af_bin  = ntile(log_MAF, 30),
    enriched = odds_ratio > 5 & padj < 0.05,
    present_afr = gvs_afr_af > 0, present_eur = gvs_eur_af > 0,
    present_eas = gvs_eas_af > 0, present_amr = gvs_amr_af > 0,
    present_sas = gvs_sas_af > 0, present_mid = gvs_mid_af > 0,
    present_oth = gvs_oth_af > 0
  ) %>%
  mutate(shared_all_continents = present_afr & present_eur & present_eas &
           present_amr & present_mid & present_sas & present_oth)

# AF bins -> representative x (geometric mean of the bin's AF range).
bin_key <- AncestrySkewAnnotated %>%
  group_by(af_bin) %>%
  summarise(min_log = min(log_MAF, na.rm = TRUE),
            max_log = max(log_MAF, na.rm = TRUE), .groups = "drop") %>%
  mutate(min_af = 10^(-max_log), max_af = 10^(-min_log), x = sqrt(min_af * max_af)) %>%
  select(af_bin, x)

all_pops   <- sort(unique(AncestrySkewAnnotated$gvs_max_subpop))
bin_totals <- AncestrySkewAnnotated %>% count(af_bin, name = "total_variants")

# Enriched-variant composition by ancestry, per AF bin.
stack_df <- AncestrySkewAnnotated %>%
  filter(enriched) %>%
  count(af_bin, gvs_max_subpop, name = "n_enriched") %>%
  complete(af_bin = bin_key$af_bin, gvs_max_subpop = all_pops, fill = list(n_enriched = 0)) %>%
  left_join(bin_totals, by = "af_bin") %>%
  mutate(total_variants = replace_na(total_variants, 0),
         prop = if_else(total_variants > 0, n_enriched / total_variants, 0)) %>%
  left_join(bin_key, by = "af_bin")

# Proportion of variants shared across all continents, per AF bin (smoothed).
line_df <- AncestrySkewAnnotated %>%
  group_by(af_bin) %>%
  summarise(prop_shared = mean(shared_all_continents, na.rm = TRUE), .groups = "drop") %>%
  complete(af_bin = bin_key$af_bin, fill = list(prop_shared = 0)) %>%
  left_join(bin_key, by = "af_bin") %>%
  arrange(x) %>%
  mutate(prop_shared_smooth = slide_dbl(prop_shared, mean, .before = 2, .after = 2, .complete = FALSE))

# Save the data for an assembly/re-plot script.
saveRDS(list(stack_df = stack_df, line_df = line_df),
        file.path(DATA_DIR, "S_ancestry_stacked.rds"))

stacked_histogram <- stack_df %>%
  mutate(pop = pretty_pop(gvs_max_subpop)) %>%
  ggplot(aes(x, prop, fill = pop)) +
    geom_area(position = "stack", alpha = 0.95) +
    geom_smooth(inherit.aes = FALSE, data = line_df, aes(x, prop_shared_smooth),
                color = "black", linewidth = 0.4, se = FALSE) +
    geom_vline(xintercept = 0.5, color = "white") +
    scale_fill_carto_d(palette = "Temps") +
    scale_x_log10(breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 5e-1),
                  labels = c("0.0001", "0.001", "0.01", "0.1", "0.5"), expand = c(0, 0)) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       expand = c(0, 0.001), breaks = c(.25, .5, .75, 1)) +
    labs(x = "AoU allele frequency",
         y = "Proportion of ancestry-enriched variants",
         fill = "Enriched\nancestry") +
    theme(axis.line = element_line(color = "black", linewidth = 0.3),
          axis.ticks = element_line(linewidth = 0.2),
          axis.title.y = element_text(margin = margin(r = 20)),
          legend.position = c(0.9, 0.5))

saveRDS(stacked_histogram, file.path(OBJ_DIR, "S_ancestry_stacked.rds"))

# Standalone supplemental figure (single panel; sized like a main-figure panel
# so the shared text sizes read correctly).
ggsave(file.path(OUTDIR, "FigS_AncestryStacked.pdf"), stacked_histogram,
       width = 13, height = 9, dpi = DPI)
showtext_opts(dpi = 100)
ggsave(file.path(OUTDIR, "FigS_AncestryStacked.png"), stacked_histogram,
       width = 13, height = 9, dpi = 100)
showtext_opts(dpi = DPI)

message("Wrote FigS_AncestryStacked.pdf / .png")
