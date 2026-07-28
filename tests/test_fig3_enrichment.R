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

cat("test_fig3_enrichment.R: all assertions passed\n")
