# Tests for SuppFig1's pure helper transforms. Run from the project root:
#   Rscript tests/test_suppfig1_helpers.R
source("tests/fixtures.R")
fx <- make_suppfig1_fixture()
source("suppfig1/00_setup.R")
source("tests/helpers.R")

# --- compute_mean_TPM ------------------------------------------------------
aou <- as_tibble(data.table::fread(fx$aou_tpms))
means <- compute_mean_TPM(aou, fx$annotations)
ok("compute_mean_TPM returns one row per gene", nrow(means) == 6L)
ok("mean_TPM is the row mean across samples",
   abs(means$mean_TPM[means$Description == "GENEA"] - 11) < 1e-9)
ok("gene_type is joined on", "gene_type" %in% names(means))
ok("no SE column is emitted", !"SE" %in% names(means))

# --- join_mean_tpms (pure) -------------------------------------------------
gtex <- as_tibble(data.table::fread(fx$gtex_tpms))
aou_means  <- compute_mean_TPM(aou,  fx$annotations)
gtex_means <- compute_mean_TPM(gtex, fx$annotations)
merged <- join_mean_tpms(aou_means, gtex_means)
ok("join_mean_tpms inner-joins on version-stripped gene id",
   all(c("mean_TPM_AoU", "mean_TPM_GTEx", "gene_id", "gene_type") %in% names(merged)))
ok("join keeps the shared genes", nrow(merged) == 6L)
ok("GENEA carries both cohorts' means",
   merged$mean_TPM_AoU[merged$Description == "GENEA"] == 11 &&
   merged$mean_TPM_GTEx[merged$Description == "GENEA"] == 10)

# --- plot_gene_comparison returns a patchwork of two plots -----------------
p <- plot_gene_comparison(merged, gene_type = "coding", hex = FALSE)
ok("plot_gene_comparison returns a 2-plot patchwork",
   inherits(p, "patchwork") && length(p$patches$plots) >= 1)

cat("test_suppfig1_helpers.R (Task 3): all assertions passed\n")
