# %% Source ancestry methods
#| Recalculate the ancestry-frequency summaries from the combined expression
#| SuSiE result table. The source result rows contain the ancestry OR and P value.
# @include preparation/Fig2/methods.R
# @include preparation/Fig2/current_methods.R
# @include notebooks/current/plots/Fig2.R

# %% Read expression results
#| Validate and hash the original table. Retain source ancestry codes.
prep_result_runtime("Fig2")
ancestry_columns <- c("pip", "gvs_max_subpop", "gvs_odds_ratio", "gvs_p_value", "gvs_all_af", "gvs_max_af",
                      "gvs_afr_af", "gvs_amr_af", "gvs_eas_af", "gvs_eur_af", "gvs_mid_af", "gvs_sas_af", "gvs_oth_af")
expression <- as_tibble(data.table::fread(input_path("susie_eqtl"), select = ancestry_columns,
                                         nThread = 4, showProgress = FALSE))
message("Read combined expression SuSiE results for S16.")

# %% Calculate ancestry-frequency summaries
#| Keep PIP >0.9 without variant deduplication. Use BH FDR over those result
#| rows. Enrichment requires OR >5 and FDR <0.05. Split -log10(MAF) into
#| 30 equal-count bins. Each x coordinate is the geometric mean of the bin's
#| AF limits. Divide each ancestry's enriched count by all rows in that bin.
#| Shared variants have nonzero AF in all seven groups. Calculate a centered
#| rolling mean using two bins on either side, with shorter windows at the ends.
#| Retain the missing-AF group in the tables. Its x coordinate is undefined,
#| so the plot layers omit that group.
bins <- fig2_ancestry_bins(expression)
stopifnot(all(dplyr::between(bins$stack_df$prop, 0, 1)))
message("Calculated 30 frequency bins, the missing-AF group, enriched counts, and shared-variant summaries.")

# %% Plot ancestry-frequency curves
#| The black LOESS curve is a visual guide through the rolling mean. It is
#| unconstrained and can exceed 100%. Display the source OTH group as REM.
figure <- fig2_plot_s16(bins)
print(figure)

# %% Save aggregate outputs
#| Save the plot and aggregate variant counts. No participant records are used.
save_plot(figure, "supplements/S16", 184, 125)
readr::write_tsv(figure$data %>% select(af_bin, gvs_max_subpop, n_enriched, total_variants, prop, x) %>%
                   mutate(gvs_max_subpop = recode(as.character(gvs_max_subpop), OTH = "REM", oth = "REM")),
                 output_file("supplements/S16_frequency_bins.tsv"))
