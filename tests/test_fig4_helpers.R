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

cat("test_fig4_helpers.R (Task 2): all assertions passed\n")
