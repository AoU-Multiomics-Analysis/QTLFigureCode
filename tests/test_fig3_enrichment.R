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
  variant     = paste0("v", 1:3),
  consequence = c("missense_variant", "intron_variant", "."),
  phylop      = c(3.5, 0.2, 1.0),
  distTSS     = c(500, 50000, -100),
  gene_type   = c("protein_coding", "lncRNA", "protein_coding"),
  max_af      = c(0.001, 0.20, 0.05),
  PLS         = c(TRUE, FALSE, TRUE),
  stringsAsFactors = FALSE))

ok("rare uses MAF < 0.01",        identical(fmx$rare, c(TRUE, FALSE, FALSE)))
ok("conserved uses phyloP > 2",   identical(fmx$conserved, c(TRUE, FALSE, FALSE)))
ok("tss_proximal uses |distTSS| < 2000 and is sign-agnostic",
   identical(fmx$tss_proximal, c(TRUE, FALSE, TRUE)))
ok("protein_coding reads gene_type",
   identical(fmx$protein_coding, c(TRUE, FALSE, TRUE)))
ok("VEP terms are expanded into the matrix", "missense_variant" %in% names(fmx))
ok("feature_names() returns only logical columns",
   all(vapply(fmx[feature_names(fmx)], is.logical, logical(1))))
ok("feature_names() excludes the outcome-ish columns",
   !any(c("variant", "consequence", "gene_type") %in% feature_names(fmx)))

cat("test_fig3_enrichment.R: all assertions passed\n")
