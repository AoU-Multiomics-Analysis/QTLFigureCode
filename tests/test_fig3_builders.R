# Tests for Fig3's derived-data builders against the mock fixture.
# Run from the project root:  Rscript tests/test_fig3_builders.R
source("tests/fixtures.R")
fx <- make_fig3_fixture()
use_fig3_fixture(fx)

source("fig3/00_setup.R")
source("tests/helpers.R")

# --- build_or_cache -------------------------------------------------------
calls <- 0
counted <- function(refresh = FALSE)
  build_or_cache("counted", refresh = refresh, { calls <<- calls + 1; calls })

ok("first call evaluates the expression", counted() == 1L)
ok("second call reads the cache",         counted() == 1L)
ok("the expression ran exactly once",     calls == 1L)
ok("refresh = TRUE rebuilds",             counted(refresh = TRUE) == 2L)
ok("cache file is written to DATA_DIR",   file.exists(file.path(DATA_DIR, "counted.rds")))

# --- susie_finemapping ----------------------------------------------------
fm <- susie_finemapping()

ok("returns all three assays plus derived tables",
   all(c("eQTL", "sQTL", "pQTL", "cs_list", "gene_map", "genes_all_three_assays") %in% names(fm)))
ok("eQTL trait ids are version-stripped",
   setequal(fm$eQTL$molecular_trait_id, c("ENSG001", "ENSG002", "ENSG003")))
ok("sQTL trait ids reduce to the gene",
   setequal(fm$sQTL$molecular_trait_id, c("ENSG001", "ENSG002")))
ok("pQTL drops the non-COMB row",          nrow(fm$pQTL) == 9L)
# Row 1 is ENSG001, whose credible set holds the 5 odd-indexed variants
# (pip 0.9, 0.7, 0.5, 0.3, 0.1).
ok("cpip sums pip within a credible set",  abs(fm$eQTL$cpip[1] - 2.5) < 1e-9)
ok("number_variants counts the credible set", fm$eQTL$number_variants[1] == 5L)
ok("cs_list labels every assay",           setequal(fm$cs_list$type, c("eQTL", "pQTL", "sQTL")))

# genes_all_three_assays is an INTERSECTION; the original misnamed it QTL_union_gene_list.
ok("genes_all_three_assays is the three-way intersection",
   setequal(fm$genes_all_three_assays, c("ENSG001", "ENSG002")))
ok("eQTL-only genes are excluded",         !"ENSG003" %in% fm$genes_all_three_assays)
ok("gene_map carries the UniProt id",      "P00001" %in% fm$gene_map$UniProt)

cat("test_fig3_builders.R: all assertions passed\n")
