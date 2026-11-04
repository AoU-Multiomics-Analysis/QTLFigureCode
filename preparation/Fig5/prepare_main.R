# Build the current five-panel Figure 5 inputs from the selected all-assay rows.
source("preparation/common.R")
source("preparation/Fig5/functions.R")
config <- prep_init()
context <- readRDS(prep_input("source", "Fig5/preparation/all_assay_context.rds"))
raw <- context$raw
reference <- context$reference
pairs <- context$pairs
frequencies <- raw %>% select(variant, all_of(fig5_af_columns)) %>% distinct()
assignments <- fig5_variant_assignments(pairs, frequencies, reference)
summary <- fig5_counts_and_classes(assignments, pairs)
funnel <- fig5_funnel(raw, reference)
effects <- fig5_effect_panels(raw, funnel, reference)
coverage_tests <- fig5_coverage_tests(summary$coverage, assignments)
primary <- bind_rows(
  fig5_subset_test(funnel$sets$rare_all, funnel$sets$pip09,
    "pip09", "rare_all", "PIP > 0.9 versus starting pool"),
  fig5_subset_test(funnel$sets$rare_all, funnel$sets$posterior_lt_neg03,
    "posterior_lt_neg03", "rare_all", "PIP > 0.9 and effect < -0.3 versus starting pool")
)
incremental <- fig5_subset_test(funnel$sets$pip09, funnel$sets$posterior_lt_neg03,
  "posterior_lt_neg03", "pip09", "Additional effect filter within high-PIP variants")

# Read original summaries only. Never compare a rerun with its own output overlay.
compare_saved <- function(observed, relative, keys, columns = names(observed)) {
  character_columns <- names(observed)[vapply(observed, is.character, logical(1))]
  column_types <- do.call(readr::cols, setNames(
    rep(list(readr::col_character()), length(character_columns)), character_columns))
  expected <- read_tsv(prep_reference("paper", relative), col_types = column_types,
    show_col_types = FALSE)
  expected <- expected %>% select(all_of(columns)) %>%
    mutate(across(where(is.factor), as.character)) %>% arrange(across(all_of(keys)))
  observed <- observed %>% select(all_of(columns)) %>%
    mutate(across(where(is.factor), as.character)) %>% arrange(across(all_of(keys)))
  if (!isTRUE(all.equal(observed, expected, check.attributes = FALSE, tolerance = 1e-10))) {
    stop("Prepared values differ from the accepted summary: ", relative, call. = FALSE)
  }
  message("Matched accepted summary: ", relative)
}
base <- "output/figure5_combined_draft"
print_base <- paste0(base, "/science_print_size")
rem_base <- paste0(print_base, "/stacked_left_arrows_REM")
funnel_base <- "output/clingen_inheritance/ar_pathogenic_unclassified_funnel"
compare_saved(summary$counts, paste0(print_base, "/data/A_counts.tsv"),
  c("allelicRequirements", "population"))
compare_saved(summary$proportions, paste0(print_base, "/data/A_clinvar_proportions.tsv"),
  c("allelicRequirements", "clinvar_class"))
compare_saved(summary$coverage, paste0(print_base, "/data/B_coverage.tsv"),
  c("allelicRequirements", "population"))
compare_saved(funnel$stages, paste0(funnel_base, "/stage_counts.tsv"), "stage")
compare_saved(effects$D, paste0(rem_base, "/data/D_display.tsv"), "variant")
compare_saved(effects$E, paste0(rem_base, "/data/E_display.tsv"), "variant")
compare_saved(coverage_tests$pooled, paste0(base,
  "/panel_B_fisher_tests/pooled_non_EUR_vs_EUR.tsv"), "allelicRequirements")
compare_saved(coverage_tests$groups, paste0(base,
  "/panel_B_fisher_tests/individual_groups_vs_EUR.tsv"), c("allelicRequirements", "comparison"))
compare_saved(coverage_tests$overall, paste0(base,
  "/panel_B_fisher_tests/unique_variants_non_EUR_vs_EUR.tsv"), "allelicRequirements")
compare_saved(coverage_tests$all, paste0(base,
  "/panel_B_fisher_tests/all_tests.tsv"), c("allelicRequirements", "comparison"))
compare_saved(primary, paste0(base, "/stage_enrichment_tests/primary_enrichment_tests.tsv"),
  "selected_stage")
compare_saved(incremental, paste0(base, "/stage_enrichment_tests/secondary_incremental_test.tsv"),
  "selected_stage")
stopifnot(n_distinct(assignments$variant) == 4479L,
  nrow(effects$D) == 14L, nrow(effects$E) == 16L,
  all(effects$D$pip > 0.9), all(effects$E$pip > 0.9))

write_paper <- function(data, path) prep_write_tsv(data, "paper", path)
write_paper(assignments, paste0(base, "/panels/AB_variant_assignments_v16.tsv"))
write_paper(frequencies %>% filter(variant %in% assignments$variant),
  paste0(base, "/panels/AB_allele_frequencies_v16.tsv"))
write_paper(summary$counts, paste0(base, "/panels/A_counts_v17_data.tsv"))
write_paper(summary$classes, paste0(base, "/panels/A_clinvar_variant_classes_v21.tsv"))
write_paper(summary$proportions, paste0(base, "/panels/A_clinvar_proportions_v21.tsv"))
write_paper(summary$coverage, paste0(base, "/panels/B_clinvar_coverage_v16_data.tsv"))
write_paper(summary$counts, paste0(print_base, "/data/A_counts.tsv"))
write_paper(summary$proportions, paste0(print_base, "/data/A_clinvar_proportions.tsv"))
write_paper(summary$coverage, paste0(print_base, "/data/B_coverage.tsv"))
write_paper(funnel$sets$rare_all,
  "output/clingen_inheritance/ar_posterior_splicing_comparison/rare_all_pairs.tsv")
write_paper(funnel$stages, paste0(funnel_base, "/stage_counts.tsv"))
write_paper(funnel$audit, paste0(funnel_base, "/gene_filter_audit.tsv"))
write_paper(funnel$shared_genes, paste0(funnel_base, "/shared_genes.tsv"))
write_paper(funnel$selected_pairs, paste0(funnel_base, "/selected_variant_gene_pairs.tsv"))
write_paper(effects$D, paste0(rem_base, "/data/D_display.tsv"))
write_paper(effects$E, paste0(rem_base, "/data/E_display.tsv"))
write_paper(effects$benign_audit,
  "output/clingen_inheritance/hexa_pmm2_galm_effects/benign_pip09_variants.tsv")
write_paper(effects$comparisons,
  "output/clingen_inheritance/ar_posterior_threshold_comparison/pathogenic_effect_comparisons_current.tsv")
iwalk(coverage_tests, function(data, name) {
  filenames <- c(pooled = "pooled_non_EUR_vs_EUR.tsv", groups = "individual_groups_vs_EUR.tsv",
    overall = "unique_variants_non_EUR_vs_EUR.tsv", all = "all_tests.tsv")
  write_paper(data, paste0(base, "/panel_B_fisher_tests/", filenames[[name]]))
})
write_paper(primary, paste0(base, "/stage_enrichment_tests/primary_enrichment_tests.tsv"))
write_paper(incremental, paste0(base, "/stage_enrichment_tests/secondary_incremental_test.tsv"))
write_paper(funnel$stages %>% slice_head(n = 3) %>%
  select(stage, variants, genes, pathogenic_variants),
  paste0(base, "/stage_enrichment_tests/verified_stage_counts.tsv"))
write_paper(imap_dfr(funnel$sets, ~ .x %>% distinct(variant, is_pathogenic) %>%
  mutate(stage = .y, .before = 1)),
  paste0(base, "/stage_enrichment_tests/variant_stage_membership.tsv"))
message("Prepared A-E and the accepted aggregate tests. Original inputs remain unchanged.")
