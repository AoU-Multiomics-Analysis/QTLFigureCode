source("preparation/common.R")
source("preparation/Fig2/methods.R")
source("preparation/Fig2/current_methods.R")
source("preparation/Fig3/result_io.R")
prep_result_runtime("Fig2")

# Read saved SuSiE results, not participant records or raw association data.
read_qtl <- function(key) as_tibble(data.table::fread(input_path(key), nThread = 4, showProgress = FALSE))
message("Read the three combined-assay SuSiE result tables")
assays <- list(Expression = read_qtl("susie_eqtl"), Splicing = read_qtl("susie_sqtl"), Proteomics = read_qtl("susie_pqtl"))

# Current A uses the separately retained frequency-bin source tables. Their
# splicing export has a historical double dot. Keep that source choice explicit.
frequency_paths <- c(Expression = "COMB.expression_SusieMerged.annotated.ancestry_skew.tsv.gz",
                     Proteomics = "COMB.proteomics_SusieMerged.annotated.ancestry_skew.tsv.gz",
                     Splicing = "COMB.splicing_SusieMerged..annotated.ancestry_skew.tsv.gz")
frequency_inputs <- map(frequency_paths, function(path) {
  alternatives <- unique(c(path, sub("SusieMerged[.][.]annotated", "SusieMerged.annotated", path)))
  source <- prep_input_any("susie", alternatives)
  readr::read_tsv(source, col_types = cols_only(variant = col_character(), molecular_trait_id = col_character(),
                                              pip = col_double(), gvs_max_af = col_double()), progress = FALSE)
})
frequency <- fig2_frequency_counts(frequency_inputs)
prep_write_tsv(frequency$counts, "paper", "output/frequency_bins/frequency_bin_counts.tsv")
prep_write_tsv(frequency$qc, "paper", "output/frequency_bins/input_qc.tsv")

# Retain the original A cache as well as the current four-bin summary.
causal <- imap_dfr(assays[c("Proteomics", "Expression", "Splicing")], ~transmute(.x, variant, pip, gvs_max_af, posterior_mean, Assay = .y)) %>%
  filter(pip > 0.9, !is.na(gvs_max_af), gvs_max_af > 0)
prep_save_rds(list(df = causal, x_center = mean(c(min(log10(causal$gvs_max_af)), log10(0.01)))), "source", "Fig2/data/A.rds")
ancestry <- fig2_ancestry_data(assays)
prep_save_rds(ancestry, "source", "Fig2/data/B.rds")
prep_save_rds(fig2_afc_data(assays$Expression), "source", "Fig2/data/C.rds")

# The historical source selects defined aFC. The November 2 audit showed that
# every defined value in these inputs already has PIP >0.9. Do not add a new
# selection rule silently. Directional effects are ordered within each gene.
cumulative <- assays$Expression %>% filter(!is.na(log2_aFC)) %>%
  make_cumulative_afc_curve_both(gene_col = "molecular_trait_id")
prep_save_rds(cumulative, "source", "Fig2/data/D.rds")
prep_write_tsv(fig2_cumulative_display(cumulative)$label_df, "paper", "output/figure2D_recalculation_20261102/prepared_display_percentages.tsv")

# Recalculate E statistics from the saved contingency counts. The earlier
# association/outlier classification is upstream of this preparation boundary.
outliers <- read_qtl("burden_enrichment")
recomputed <- fig2_outlier_statistics(outliers)
for (column in c("odds_ratio", "se_log_or", "ci_low", "ci_high", "p_value")) {
  stopifnot(isTRUE(all.equal(outliers[[column]], recomputed[[column]], tolerance = 1e-10)))
}
prep_save_rds(recomputed, "source", "Fig2/data/E.rds")
prep_save_rds(fig2_median_data(read_qtl("median_ko")), "source", "Fig2/data/F.rds")

# S16: original 30 equal-count -log10(MAF) bins, BH test universe, geometric
# bin coordinates, seven-population presence flag, and centered rolling mean.
bins <- fig2_ancestry_bins(assays$Expression)
prep_save_rds(bins, "source", "Fig2/data/S_ancestry_stacked.rds")
prep_write_tsv(bins$stack_df, "paper", "output/figure2_ancestry_bins/stacked_counts.tsv")
prep_write_tsv(bins$line_df, "paper", "output/figure2_ancestry_bins/shared_curve.tsv")
message("Figure 2 A–F and S16 plotting data prepared; upstream QTL models were not rerun")
