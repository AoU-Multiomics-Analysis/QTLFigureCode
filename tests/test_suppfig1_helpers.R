# Tests for SuppFig1's pure helper transforms. Run from the project root:
#   Rscript tests/test_suppfig1_helpers.R
source("tests/fixtures.R")
fx <- make_suppfig1_fixture()
source("suppfig1/00_setup.R")
source("tests/helpers.R")

# --- compute_mean_TPM ------------------------------------------------------
# The AoU GCT has 7 genes: 6 shared with GTEx plus 1 AoU-only (see fixtures.R).
aou <- as_tibble(data.table::fread(fx$aou_tpms))
means <- compute_mean_TPM(aou, fx$annotations)
ok("compute_mean_TPM returns one row per gene", nrow(means) == 7L)
ok("mean_TPM is the row mean across samples",
   abs(means$mean_TPM[means$Description == "GENEA"] - 11) < 1e-9)
ok("gene_type is joined on", "gene_type" %in% names(means))
ok("no SE column is emitted", !"SE" %in% names(means))

# --- join_mean_tpms (pure) -------------------------------------------------
# The GTEx GCT only has the 6 shared genes (no AoU-only gene), so the inner
# join in join_mean_tpms must drop the AoU-only gene from the merged result.
gtex <- as_tibble(data.table::fread(fx$gtex_tpms))
aou_means  <- compute_mean_TPM(aou,  fx$annotations)
gtex_means <- compute_mean_TPM(gtex, fx$annotations)
merged <- join_mean_tpms(aou_means, gtex_means)
shared_genes <- c("GENEA", "TFGENE", "LNC1", "SNR1", "HBB", "RRNA1")  # 6 shared genes
ok("join_mean_tpms inner-joins on version-stripped gene id",
   all(c("mean_TPM_AoU", "mean_TPM_GTEx", "gene_id", "gene_type") %in% names(merged)))
ok("join keeps exactly the 6 shared genes (AoU-only gene dropped by inner join)",
   nrow(merged) == 6L)
ok("the AoU-only gene is absent from the merged result",
   !("AOUONLY7" %in% merged$Description))
ok("all shared genes are present in the merged result",
   setequal(merged$Description, shared_genes))
ok("GENEA carries both cohorts' means",
   merged$mean_TPM_AoU[merged$Description == "GENEA"] == 11 &&
   merged$mean_TPM_GTEx[merged$Description == "GENEA"] == 10)

# --- plot_gene_comparison returns a patchwork of two plots -----------------
p <- plot_gene_comparison(merged, gene_type = "coding", hex = FALSE)
ok("plot_gene_comparison returns a 2-plot patchwork",
   inherits(p, "patchwork") && length(p$patches$plots) >= 1)

# --- cohort_specific_genes (pure) ------------------------------------------
spec_mock <- tibble::tibble(
  Description = c("A", "B", "C", "D"),
  gene_type = "protein_coding",
  mean_TPM_AoU  = c(5,   0.2, 5,   0.2),
  mean_TPM_GTEx = c(0.2, 5,   5,   0.2))
ok("AoU-specific = coding, AoU>1 & GTEx<1",
   identical(cohort_specific_genes(spec_mock, "AoU"), "A"))
ok("GTEx-specific = coding, AoU<1 & GTEx>1",
   identical(cohort_specific_genes(spec_mock, "GTEx"), "B"))

cat("test_suppfig1_helpers.R (Task 3): all assertions passed\n")
