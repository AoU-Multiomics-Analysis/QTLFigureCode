# %% Source result methods
#| Start from the combined-assay SuSiE result tables. These cells show the
#| frequency filters, ancestry tests, annotation selection, aFC models, and
#| contingency-table tests used for Figure 2. The QTL models are upstream inputs.
# @include preparation/Fig2/methods.R
# @include preparation/Fig2/current_methods.R
# @include notebooks/current/plots/Fig2.R

# %% Read the source result tables
#| Read expression, splicing, and protein QTL results. Keep each result row
#| until the panel-specific selection below. Input paths are checked and hashed.
prep_result_runtime("Fig2")
read_qtl <- function(key, columns = NULL) {
  as_tibble(data.table::fread(input_path(key), select = columns, nThread = 4, showProgress = FALSE))
}
ancestry_columns <- c("variant", "pip", "gvs_max_subpop", "gvs_max_af", "gvs_odds_ratio", "gvs_p_value")
expression_columns <- c(ancestry_columns, "log2_aFC", "posterior_mean", "molecular_trait_id", "consequence",
                        "FANTOM5", "pELS", "dELS", "PLS", "CTCF.bound", "DNase.H3K4me3", "CTCF.only")
assays <- list(Expression = read_qtl("susie_eqtl", expression_columns),
               Splicing = read_qtl("susie_sqtl", ancestry_columns),
               Proteomics = read_qtl("susie_pqtl", ancestry_columns))
message("Read three combined-assay SuSiE result tables.")

# %% Four frequency bins
#| A counts high-PIP variant-trait pairs. The denominator is specific to each
#| assay. Retain the current <1%, 1–<5%, 5–<10%, and 10–50% bins.
#| A uses the separately retained source tables. The splicing file has a
#| historical double dot. The alternative spelling is checked when needed.
frequency_paths <- c(Expression = "COMB.expression_SusieMerged.annotated.ancestry_skew.tsv.gz",
                     Proteomics = "COMB.proteomics_SusieMerged.annotated.ancestry_skew.tsv.gz",
                     Splicing = "COMB.splicing_SusieMerged..annotated.ancestry_skew.tsv.gz")
frequency_inputs <- map(frequency_paths, function(path) {
  candidates <- unique(c(path, sub("SusieMerged[.][.]annotated", "SusieMerged.annotated", path)))
  readr::read_tsv(prep_input_any("susie", candidates),
                  col_types = cols_only(variant = col_character(), molecular_trait_id = col_character(),
                                        pip = col_double(), gvs_max_af = col_double()), progress = FALSE)
})
frequency <- fig2_frequency_counts(frequency_inputs)
plots <- list(A = fig2_plot_a(frequency$counts))
readr::write_tsv(frequency$counts, output_file("Fig2/frequency_bin_counts.tsv"))
readr::write_tsv(frequency$qc, output_file("Fig2/frequency_input_qc.tsv"))
save_plot(plots$A, "Fig2/panel_A", 150, 105)
print(plots$A)

# %% Ancestry enrichment
#| B keeps PIP >0.9 and one row per variant across assays. Recalculate BH FDR
#| after that selection. An enriched variant has OR >5 and FDR <0.05.
#| OR and P values come from the source ancestry-enrichment tests. Rare means
#| population-max AF <1%. OTH is displayed as REM.
ancestry <- fig2_ancestry_data(assays)
plots$B <- fig2_plot_b(ancestry)
save_plot(plots$B, "Fig2/panel_B", 150, 105)
print(plots$B)

# %% Allelic fold change
#| C selects the highest-priority retained consequence for each high-PIP
#| expression variant. Fit |log2 aFC| ~ rare within each consequence group and
#| apply BH FDR across those estimates. The plot shows medians and quantiles.
#| The 0–2 view clips the distributions but does not remove source values.
afc <- fig2_afc_data(assays$Expression)
plots$C <- fig2_plot_c(afc)
c_artwork <- render_panel_C(plots$C)
save_plot(c_artwork, "Fig2/panel_C", 150, 105)
grid::grid.newpage()
grid::grid.draw(c_artwork)

# %% Cumulative effects
#| D selects defined aFC values, as in the source calculation. In the current
#| input, every defined value already has PIP >0.9. Order decreasing and
#| increasing effects within each gene, then sum log2 aFC in each direction.
#| Expected heterozygous expression is (1 + 2^cumulative_log2_aFC) / 2.
#| Retain values below 400% before pooling steps into 4+. The label denominator
#| is the displayed gene-step observations. The thresholds are 60% and 150%.
cumulative <- assays$Expression %>% filter(!is.na(log2_aFC)) %>%
  make_cumulative_afc_curve_both(gene_col = "molecular_trait_id")
cumulative_display <- fig2_cumulative_display(cumulative)
plots$D <- fig2_plot_d(cumulative_display)
readr::write_tsv(cumulative_display$label_df, output_file("Fig2/cumulative_display_percentages.tsv"))
save_plot(plots$D, "Fig2/panel_D", 150, 105)
print(plots$D)

# %% Outlier contingency tests
#| E starts from the saved aggregate contingency counts. The association and
#| outlier classifications are upstream inputs. Recalculate each OR with 0.5
#| added to all four cells, its Wald interval, and a two-sided Fisher P value
#| from the original integer counts. Omit the central reference bin from E.
outlier_counts <- read_qtl("burden_enrichment")
outliers <- fig2_outlier_statistics(outlier_counts)
for (column in c("odds_ratio", "se_log_or", "ci_low", "ci_high", "p_value")) {
  stopifnot(isTRUE(all.equal(outlier_counts[[column]], outliers[[column]], tolerance = 1e-10)))
}
message("Recalculated outlier odds ratios, confidence intervals, and Fisher tests.")

# %% Gene burden
#| F starts from aggregate per-bin medians and quartiles. Keep protein-coding
#| genes, correct the recorded central-bin label, and derive each lower bound.
#| Participant records are not used. The tail annotations sum bin medians.
median <- fig2_median_data(read_qtl("median_ko"))
plots$E <- fig2_plot_e(outliers, median)
plots$F <- fig2_plot_f(median)
bin_lower <- sort(unique(median$lower))
bin_transform <- fig2_bin_transform(bin_lower)
e_build <- ggplot_build(plots$E)
f_build <- ggplot_build(plots$F)
e_range <- e_build$layout$panel_params[[1]]$x.range
f_range <- f_build$layout$panel_params[[1]]$x.range
alignment <- tibble(lower_bound = sort(unique(plots$E$data$focal_lower_bound))) %>%
  mutate(bin_position = match(lower_bound, bin_lower),
         E_fraction = (bin_transform$transform(lower_bound) - e_range[[1]]) / diff(e_range),
         F_fraction = (bin_position - f_range[[1]]) / diff(f_range),
         aligned = abs(E_fraction - F_fraction) < 1e-10)
stopifnot(all(alignment$aligned))
readr::write_tsv(alignment, output_file("Fig2/EF_alignment_check.tsv"))
ef <- (plots$E / plots$F) + plot_layout(heights = c(24, 33))
save_plot(ef, "Fig2/panels_EF", 183, 90)
print(ef)

# %% Complete Figure 2
#| Assemble the six freshly constructed panels at the current print width.
layout_plots <- plots
layout_plots$C <- wrap_elements(full = c_artwork)
fig2 <- wrap_plots(layout_plots, design = "AB\nCD\nEE\nFF", widths = c(81, 98),
                   heights = c(48, 62, 24, 33)) + plot_annotation(tag_levels = "A")
save_plot(fig2, "Fig2/Figure2_rebuilt", 183, 188)
print(fig2)

# %% Reference figure
#| Complete figure layout.
show_reference(2)
