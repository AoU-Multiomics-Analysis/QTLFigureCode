# Tests for Fig4's pure helper transforms. Run from the project root:
#   Rscript tests/test_fig4_helpers.R
source("tests/fixtures.R")
fx <- make_fig4_fixture(); use_fig4_fixture(fx)
source("fig4/00_setup.R")
source("tests/helpers.R")

# --- build_plot_data -------------------------------------------------------
gene <- as_tibble(data.table::fread(input_path("coloc_gene")))
rate <- as_tibble(data.table::fread(input_path("coloc_rate")))

combined <- build_plot_data(gene, rate, min_coloc_genes = 1, rate_type = "lenient")
ok("build_plot_data returns one row per trait", nrow(combined) == 2L)
ok("carries the selected union rate",
   all(c("pct_coloc_union", "count") %in% names(combined)))
ok("count is the number of colocalizing genes per trait",
   all(combined$count == 8L))
ok("lenient rate is selected", combined$pct_coloc_union[combined$trait == "eosinophil"] == 45)

expect_error("an unknown rate_type errors",
             build_plot_data(gene, rate, rate_type = "nope"))
expect_error("a non-positive min_coloc_genes errors",
             build_plot_data(gene, rate, min_coloc_genes = 0))

# min_coloc_genes trims traits below the threshold
trimmed <- build_plot_data(gene, rate, min_coloc_genes = 9, rate_type = "lenient")
ok("min_coloc_genes trims traits with too few genes", nrow(trimmed) == 0L)

# as_logical_flag coerces common encodings
ok("as_logical_flag maps yes/1/true",
   identical(as_logical_flag(c("yes", "1", "true", "no")), c(TRUE, TRUE, TRUE, FALSE)))

# --- derive_gene_sets ------------------------------------------------------
clingen_mock <- tibble::tibble(
  targetId          = c("ENSG1", "ENSG1", "ENSG2", "ENSG3"),
  targetFromSourceId = c("S1", "S1", "S2", "S3"),
  studyId           = c("st1", "st1", "st2", "st3"),
  confidence        = c("Definitive", "Definitive", "Strong", "Limited"),
  allelicRequirements = list("AR", "AR", "AD", "AD"),
  diseaseFromSource = c("hereditary anemia", "hereditary anemia",
                        "dilated cardiomyopathy", "some limited disease"))
drug_mock <- tibble::tibble(
  gene_name = c("geneA", "geneB", "geneB", ""),
  approved  = c(TRUE, FALSE, FALSE, TRUE))

gs <- derive_gene_sets(clingen_mock, drug_mock)
ok("mendelian keeps only Definitive/Strong genes",
   setequal(gs$mendelian, c("ENSG1", "ENSG2")))
ok("Limited-confidence genes are excluded", !"ENSG3" %in% gs$mendelian)
ok("heme-immune picks the anemia gene", "ENSG1" %in% gs$heme_immune)
ok("heme-immune excludes the cardiomyopathy gene", !"ENSG2" %in% gs$heme_immune)
ok("druggable upper-cases symbols and drops blanks",
   setequal(gs$druggable, c("GENEA", "GENEB")))
ok("approved keeps only approved==TRUE non-blank symbols",
   setequal(gs$approved, "GENEA"))

# --- calculate_qq_data (pure; no MatchIt needed) ---------------------------
matched_mock <- tibble::tibble(
  group = rep(c("Colocalized genes", "Matched non-colocalized genes"), each = 4),
  coloc_trait = "eosinophil",
  common_pvalue = c(1e-8, 1e-4, 1e-2, 0.5, 0.9, 0.5, 0.2, 0.05),
  p_rare_cond   = c(1e-20, 1e-3, 0.02, 0.7, 0.8, 0.6, 0.3, 0.1))
qq <- calculate_qq_data(matched_mock, traits = "eosinophil")
ok("calculate_qq_data returns both tests", setequal(unique(qq$test), c("common", "rare")))
ok("expected/observed columns present",
   all(c("expected", "observed", "observed_display", "capped") %in% names(qq)))
ok("common p-values are capped at 40",  all(qq$observed_display[qq$test == "common"] <= 40))
ok("rare p-values are capped at 12.5",  all(qq$observed_display[qq$test == "rare"]   <= 12.5))
ok("the 1e-20 rare hit is flagged capped", any(qq$capped[qq$test == "rare"]))

# --- prepare_qq_match_source (pure filtering/slicing) ----------------------
cond_mock <- tibble::tibble(
  status = "ok",
  gene = c("g1", "g1", "g2"),
  coloc_trait = "eosinophil",
  molecular_trait_id = c("m1", "m2", "m3"),
  any_coloc = c(TRUE, TRUE, FALSE),
  z_rare_cond = c(2, 2, 1), rare_zscore = c(2, 2, 1), Vcor = c(0.1, 0.1, 0.1),
  GLCP = c(0.9, 0.8, NA), p_rare_cond = c(0.01, 0.2, 0.5),
  mean_TPM = c(10, 10, 5), shet = c(0.01, 0.01, 0.02), n_rare = c(3, 3, 1))
src <- prepare_qq_match_source(cond_mock)
ok("prepare_qq_match_source keeps one row per (gene, coloc_trait)", nrow(src) == 2L)
ok("it adds the modelling covariates",
   all(c("log_tpm", "log_shet", "log_n_rare", "abs_vcor", "coloc") %in% names(src)))

# --- match_colocalized_pairs (only if MatchIt is installed) ----------------
if (requireNamespace("MatchIt", quietly = TRUE)) {
  big <- dplyr::bind_rows(lapply(1:40, function(i) tibble::tibble(
    status = "ok", gene = paste0("g", i), coloc_trait = "eosinophil",
    molecular_trait_id = paste0("m", i),
    any_coloc = i %% 2 == 0,
    z_rare_cond = rnorm(1, 2), rare_zscore = rnorm(1, 2), Vcor = 0.1,
    GLCP = 0.9, p_rare_cond = runif(1),
    mean_TPM = runif(1, 1, 50), shet = runif(1, 0.001, 0.1), n_rare = sample(1:10, 1))))
  set.seed(1)
  mr <- match_colocalized_pairs(prepare_qq_match_source(big))
  ok("match returns a fit and a matched data frame",
     all(c("fit", "data") %in% names(mr)) && nrow(mr$data) > 0)
  ok("matched data has the two comparison groups",
     setequal(levels(mr$data$group),
              c("Matched non-colocalized genes", "Colocalized genes")))
} else {
  cat("  (skipped match_colocalized_pairs: MatchIt not installed)\n")
}

cat("test_fig4_helpers.R (Task 2): all assertions passed\n")
