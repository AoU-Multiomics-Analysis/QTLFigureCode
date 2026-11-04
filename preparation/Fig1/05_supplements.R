# Prepare the current S12 and S13 summaries and display-object input contracts.
# Fine-mapping fits are inputs. Analytic gene-cluster SEs use the archived formula.
source("preparation/common.R")
config <- prep_init()
relative_base <- "output/eas_fine_mapping_2026-11-02"
walk(c("credible_sets.tsv.gz", "lead_variant_pips.tsv.gz"),
  ~prep_input("paper", file.path(relative_base, "data", .x)))
source("preparation/Fig1/build_sizes_display.R")
source("preparation/Fig1/build_lead_display.R")
copies <- list(
  list(source = file.path(relative_base, "credible_set_sizes"),
       target = "Supplemental Figures/01_Supplemental_drafts/Figure1_population_credible_set_sizes/v5_with_EAS",
       files = c("figure_display.rds", "population_credible_set_summary.tsv", "credible_sets.tsv.gz",
                 "UKBB_reference_marker.tsv", "AFR_mean_size_differences.tsv", "gene_cluster_summaries.tsv.gz")),
  list(source = file.path(relative_base, "lead_pip"),
       target = "Supplemental Figures/01_Supplemental_drafts/Figure1_population_lead_PIP/v2_with_EAS",
       files = c("figure_display.rds", "lead_variant_pips.tsv.gz", "lead_pip_summary.tsv",
                 "lead_pip_density_curves.tsv.gz", "density_area_checks.tsv")))
walk(copies, function(entry) walk(entry$files, function(name) {
  from <- file.path(config$output_root, "paper", entry$source, name)
  stopifnot(file.copy(from, prep_output("paper", file.path(entry$target, name)), overwrite = TRUE))
}))
# Compare the actual scientific summaries with the fixed submission values.
compare_numeric <- function(folder, filename, keys) {
  generated <- read_tsv(file.path(config$output_root, "paper", folder, filename), show_col_types = FALSE)
  saved <- read_tsv(prep_reference("paper", file.path(folder, filename)), show_col_types = FALSE)
  check <- generated %>% arrange(across(all_of(keys)))
  expected <- saved %>% arrange(across(all_of(keys)))
  stopifnot(isTRUE(all.equal(check, expected, tolerance = 1e-10, check.attributes = FALSE)))
}
compare_numeric(copies[[1]]$target, "AFR_mean_size_differences.tsv", c("assay", "population"))
compare_numeric(copies[[1]]$target, "population_credible_set_summary.tsv", c("assay", "population"))
compare_numeric(copies[[2]]$target, "lead_pip_summary.tsv", c("assay", "population"))
compare_numeric(copies[[2]]$target, "lead_pip_density_curves.tsv.gz", c("assay", "population", "pip"))
message("S12/S13 preparation verified: all size summaries, analytic contrasts, PIP summaries and density values match.")
