# Fig3 panel E — colocalization feature enrichment.
# Run from the project root:  Rscript tests/test_fig3_enrichment.R
source("tests/fixtures.R")
source("tests/helpers.R")

fx <- make_fig3_fixture(file.path(tempdir(), "fig3_enrich"))
use_fig3_fixture(fx)
source("fig3/00_setup.R")

fm  <- susie_finemapping(refresh = TRUE)
map <- .variant_feature_map(fm)

ok("one row per variant", !any(duplicated(map$variant)))
ok("carries the cCRE flags", all(FEATURE_CCRE_COLS %in% names(map)))
ok("carries the annotation columns",
   all(c("consequence", "phylop", "distTSS", "gene_type") %in% names(map)))
ok("carries a per-variant median credible-set size", "cs_size_med" %in% names(map))
ok("cs_size_med is numeric and positive",
   is.numeric(map$cs_size_med) && all(map$cs_size_med > 0, na.rm = TRUE))
ok("allele frequencies survive the widening",
   all(AF_POP_COLS %in% names(map)))

# --- feature matrix ---------------------------------------------------------
fmx <- build_feature_matrix(data.frame(
  variant     = paste0("v", 1:4),
  consequence = c("missense_variant", "intron_variant", ".", "."),
  phylop      = c(3.5, 0.2, 1.0, 0.5),
  # Row 3 (-100) and row 4 (-50000) both cover the sign-agnostic claim: a
  # small-magnitude negative offset should still count as proximal, and a
  # large-magnitude negative offset should not. Row 3 alone cannot falsify an
  # implementation missing abs() -- (-100 < 2000) is TRUE either way -- but
  # row 4's (-50000 < 2000) would wrongly be TRUE without abs().
  distTSS     = c(500, 50000, -100, -50000),
  gene_type   = c("protein_coding", "lncRNA", "protein_coding", "protein_coding"),
  max_af      = c(0.001, 0.20, 0.05, 0.30),
  PLS         = c(TRUE, FALSE, TRUE, FALSE),
  stringsAsFactors = FALSE))

ok("rare uses MAF < 0.01",        identical(fmx$rare, c(TRUE, FALSE, FALSE, FALSE)))
ok("conserved uses phyloP > 2",   identical(fmx$conserved, c(TRUE, FALSE, FALSE, FALSE)))
ok("tss_proximal uses |distTSS| < 2000 and is sign-agnostic",
   identical(fmx$tss_proximal, c(TRUE, FALSE, TRUE, FALSE)))
ok("protein_coding reads gene_type",
   identical(fmx$protein_coding, c(TRUE, FALSE, TRUE, TRUE)))
ok("VEP terms are expanded into the matrix", "missense_variant" %in% names(fmx))
ok("feature_names() returns only logical columns",
   all(vapply(fmx[feature_names(fmx)], is.logical, logical(1))))
ok("feature_names() excludes the outcome-ish columns",
   !any(c("variant", "consequence", "gene_type") %in% feature_names(fmx)))

cat("test_fig3_enrichment.R: all assertions passed\n")
