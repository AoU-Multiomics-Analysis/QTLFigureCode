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

# --- .label_by_pair: FDR must be computed PER PAIR, not pooled -------------
# One pair ("strong") holds many near-certain probabilities; the other
# ("weak") holds only two mediocre ones. Under correct per-pair FDR control
# neither of the weak pair's rows (0.6, 0.5) clears a 5% threshold on its own,
# so the weak pair must come back entirely unlabelled. Pooling would let the
# strong pair's 18 near-zero local FDRs drag the pooled running mean far
# enough under 5% to also pass the weak pair's two rows -- i.e. pooling would
# have labelled both of them TRUE. The two asserts below distinguish exactly
# this: on the current mock fixture used above, per-pair and pooled happen to
# agree, so only a frame shaped like this one would catch a regression to
# pooling.
fdr_test <- tibble::tibble(
  type = c(rep("strong", 18), rep("weak", 2)),
  GLCP = c(rep(0.999, 18), 0.6, 0.5)
)
per_pair <- .label_by_pair(fdr_test, "GLCP", "fdr test")
ok("per-pair: the strong pair is fully labelled",
   all(per_pair$Coloc[per_pair$type == "strong"]))
ok("per-pair: the weak pair is left entirely unlabelled",
   !any(per_pair$Coloc[per_pair$type == "weak"]))

# Sanity check that this fixture actually distinguishes the two behaviours:
# pooling (no split by `type`) mislabels the weak pair, confirming the asserts
# above are not vacuously true regardless of per-pair vs. pooled behaviour.
pooled <- add_fdr_threshold_label(fdr_test, "GLCP", quiet = TRUE)
ok("pooling on the same frame would have wrongly labelled the weak pair",
   any(pooled$Coloc[pooled$type == "weak"]))

# --- coloc_variant_level: duplicated annotation SNP key must error --------
# The annotation join assumes SNP is unique in the (unvalidated) AoU table.
# A duplicated key must fail the build loudly (relationship = "many-to-one")
# rather than silently multiplying rows that then get cached to disk.
fx_dup_ann <- make_fig3_fixture(file.path(tempdir(), "fig3_fixture_dup_ann"))
ann_path <- file.path(fx_dup_ann$root,
                      "annotated/AoU.COMB.AC10.AN95.biallelic.v9.annotations.tsv.bgz")
ann <- read.delim(ann_path)
ann <- rbind(ann, ann[1, ])   # duplicate one SNP's row
write.table(ann, ann_path, sep = "\t", row.names = FALSE, quote = FALSE)
# INPUT/DATA_DIR/MANIFEST/OUTDIR are plain R globals, assigned once from
# Sys.getenv() when fig3/00_setup.R is sourced -- setting the env vars alone
# does not retroactively change them, so 00_setup.R must be re-sourced
# (cheap and idempotent, per its header) to pick up the new fixture.
use_fig3_fixture(fx_dup_ann)
source("fig3/00_setup.R")
expect_error("a duplicated annotation SNP key errors the variant-level join",
             coloc_variant_level(refresh = TRUE))

# --- coloc_cs_level: a Signal with more than one '(@)' must error ---------
# FastENLOC's documented Signal format is exactly two sides,
# "gene:cs1(@)gene:cs2"; a multi-credible-set side uses "_", not another
# "(@)", so two or more "(@)" indicates a malformed row that must not be
# silently folded by extra = "merge".
fx_dup_sig <- make_fig3_fixture(file.path(tempdir(), "fig3_fixture_dup_sig"))
sig_path <- file.path(fx_dup_sig$root, "enloc/eQTL_pQTL.enloc.sig.out")
sig_lines <- readLines(sig_path)
sig_lines[2] <- sub("ENSG001:P00001_L1", "ENSG001:P00001_L1(@)EXTRA", sig_lines[2])
writeLines(sig_lines, sig_path)
use_fig3_fixture(fx_dup_sig)
source("fig3/00_setup.R")
expect_error("a Signal with two '(@)' separators errors the cs-level split",
             coloc_cs_level(refresh = TRUE))

# Restore the original fixture so any later-added tests in this file (or a
# re-run) see the well-formed inputs, not the corrupted ones above.
use_fig3_fixture(fx)
source("fig3/00_setup.R")

cat("test_fig3_builders.R: all assertions passed\n")
