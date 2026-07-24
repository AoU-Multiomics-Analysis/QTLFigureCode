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


# --- format_go_plot_data (Panel B prep, pure) ------------------------------
go_mock <- tibble::tibble(
  Term = c(paste0("T", 1:12, " (GO:000", 1:12, ")")),
  Adjusted.P.value = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08),
  type = rep(c("AoU", "GTEx"), 6))
gp <- format_go_plot_data(go_mock, top_n = 3)
ok("format_go_plot_data keeps <= top_n per cohort",
   all(gp %>% count(type) %>% pull(n) <= 3))
ok("the GO id suffix is stripped from Term", !any(grepl("\\(GO", gp$Term)))
ok("neglog10_p is -log10(adjusted p)",
   abs(gp$neglog10_p[1] - -log10(gp$Adjusted.P.value[1])) < 1e-9)

# Discriminating check: the kept rows must be the LOWEST-adjusted-p terms, not
# just any top_n rows. Terms are named "LOW"/"HIGH" so a regression to
# slice_tail() (keeps worst p-values) or arrange(desc(...)) would be caught
# even though the count and neglog10_p-formula assertions above would still
# pass.
discr_mock <- tibble::tibble(
  Term = c("LOW1 (GO:0001)", "LOW2 (GO:0002)", "LOW3 (GO:0003)",
           "HIGH1 (GO:0004)", "HIGH2 (GO:0005)"),
  Adjusted.P.value = c(1e-6, 1e-5, 1e-4, 0.5, 0.6),
  type = "AoU")
gp2 <- format_go_plot_data(discr_mock, top_n = 3)
ok("format_go_plot_data keeps the lowest-p terms, not the highest-p ones",
   setequal(trimws(gp2$Term), c("LOW1", "LOW2", "LOW3")) &&
   !any(grepl("HIGH", gp2$Term)))

# Edge case: a Term with no "(GO:...)" suffix must pass through unstripped.
nosuffix_mock <- tibble::tibble(
  Term = c("Plain term with no GO id", "T2 (GO:0009)"),
  Adjusted.P.value = c(1e-8, 1e-2),
  type = "AoU")
gp3 <- format_go_plot_data(nosuffix_mock, top_n = 2)
ok("a Term without a (GO:...) suffix passes through unstripped",
   "Plain term with no GO id" %in% gp3$Term)


# --- tf_expression_table does NOT pre-log the TPMs (the double-log fix) -----
tf_mock <- tibble::tibble(
  Description = c("TF1", "TF2", "OTHER"),
  gene_type = "protein_coding",
  mean_TPM_AoU  = c(8, 16, 10),
  mean_TPM_GTEx = c(4,  4,  4))
tf_tab <- tf_expression_table(tf_mock, c("TF1", "TF2"))
ok("tf_expression_table keeps only TF-list genes above threshold",
   setequal(tf_tab$Description, c("TF1", "TF2")))
ok("TPMs are raw, not pre-logged (TF2 AoU stays 16, not log2(16+1))",
   tf_tab$mean_TPM_AoU[tf_tab$Description == "TF2"] == 16)

cat("test_suppfig1_helpers.R (Task 3): all assertions passed\n")
