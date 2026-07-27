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

# --- coloc_variant_level: conflicting allele frequencies must error --------
# Allele frequencies now come from the SuSiE tables rather than a standalone
# annotation table, so the same variant appears in up to three of them. They
# must agree: a disagreement would otherwise give one variant two AF rows,
# fanning out the variant-level join and caching inflated rows to disk.
fx_dup_ann <- make_fig3_fixture(file.path(tempdir(), "fig3_fixture_dup_ann"))
sq_path <- file.path(fx_dup_ann$root,
                     "finemapping/coloc/sQTL_susie.COMB.annotated.cleaned.tsv.gz")
sq <- read.delim(sq_path)
sq$gvs_afr_af[1] <- sq$gvs_afr_af[1] + 0.123   # disagree with eQTL/pQTL
write.table(sq, sq_path, sep = "\t", row.names = FALSE, quote = FALSE)
# INPUT/DATA_DIR/MANIFEST/OUTDIR are plain R globals, assigned once from
# Sys.getenv() when fig3/00_setup.R is sourced -- setting the env vars alone
# does not retroactively change them, so 00_setup.R must be re-sourced
# (cheap and idempotent, per its header) to pick up the new fixture.
use_fig3_fixture(fx_dup_ann)
source("fig3/00_setup.R")
expect_error("conflicting per-assay allele frequencies error the variant-level join",
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

# --- coloc_model_input ----------------------------------------------------
mi <- coloc_model_input()

ok("returns one element per assay pair",
   setequal(names(mi), c("eQTL-sQTL", "eQTL-pQTL", "sQTL-pQTL")))

ep <- mi[["eQTL-pQTL"]]
ok("adds max_af / rare / log_maf",
   all(c("max_af", "rare", "log_maf") %in% names(ep)))
ok("adds both joined posterior columns",
   all(c("posterior_mean.x", "posterior_mean.y") %in% names(ep)))
ok("keeps only colocalized rows",     all(ep$Coloc))
ok("log_maf is a plain numeric vector, not a matrix", is.null(dim(ep$log_maf)))

# chr1_100 has gvs_eur_af = 0.90, which must fold to a MAF of 0.10 before the
# per-variant maximum is taken.
chr1 <- ep[ep$SNP == "chr1_100_A_G", ]
ok("allele frequencies fold to MAF before pmax", abs(chr1$max_af[1] - 0.30) < 1e-9)
ok("common variants are not flagged rare",       isFALSE(chr1$rare[1]))

sp <- mi[["sQTL-pQTL"]]
ok("rare variants are flagged",
   nrow(sp) == 0L || all(sp$max_af[sp$rare] < 0.01))

# The above passes vacuously if `rare` is always FALSE (all() of an empty
# vector is TRUE). Pin the boundary from both sides: some rows must actually
# be flagged, the flagged ones must be below the threshold, AND the
# unflagged ones must be at or above it.
ok("some variants are actually flagged rare", any(sp$rare))
ok("unflagged variants are at or above the rare threshold",
   all(sp$max_af[!sp$rare] >= 0.01))

# --- posterior joins must not fan out row count ----------------------------
# The "adds both joined posterior columns" check above only tests that the
# column names exist, which would not have caught a fan-out join silently
# duplicating rows. Assert the join preserves row count against the
# pre-join, colocalized frame it was built from.
pre_join_rows <- coloc_variant_level() %>%
  dplyr::filter(type == "eQTL-pQTL", Coloc == TRUE) %>%
  nrow()
ok("eQTL-pQTL posterior joins preserve row count",
   nrow(ep) == pre_join_rows)

# --- coloc_model_input: duplicate-key attribution must name the CORRECT table
# pair_input() chains two guarded joins (join_posteriors). If they were
# chained through %>%, the outer join's `d` argument would be a lazy promise
# wrapping the inner join_posteriors() call; forcing it inside the outer
# tryCatch's monitored left_join() would let the outer handler catch the
# inner join's error and mis-attribute it to the OTHER table's name. For each
# of the three fine-mapping tables in turn, inject a duplicate (variant,
# gene) row and confirm the error names THAT table -- and only that table.
dup_finemap_row <- function(rel_path) {
  fx_dup <- make_fig3_fixture(file.path(tempdir(), paste0("fig3_fixture_dup_", basename(rel_path))))
  p <- file.path(fx_dup$root, rel_path)
  d <- read.delim(p)
  d <- rbind(d, d[1, ])   # duplicate row 1's (variant, gene) key
  write.table(d, p, sep = "\t", row.names = FALSE, quote = FALSE)
  use_fig3_fixture(fx_dup)
  source("fig3/00_setup.R")
  tryCatch({ coloc_model_input(refresh = TRUE); NULL },
           error = conditionMessage)
}

msg_eqtl <- dup_finemap_row("finemapping/coloc/eQTL_susie.COMB.annotated.cleaned.tsv.gz")
ok("eQTL duplicate is reported as an eQTL problem",
   grepl("^eQTL fine-mapping table has duplicate", msg_eqtl))
ok("eQTL duplicate is NOT blamed on sQTL or pQTL",
   !grepl("sQTL fine-mapping table", msg_eqtl) && !grepl("pQTL fine-mapping table", msg_eqtl))
ok("eQTL duplicate error shows an offending (variant, gene) pair",
   grepl("variant = '.*', gene = '.*'", msg_eqtl))
ok("a failed coloc_model_input build (eQTL) leaves no cache file behind",
   !file.exists(file.path(DATA_DIR, "coloc_model_input.rds")))

msg_sqtl <- dup_finemap_row("finemapping/coloc/sQTL_susie.COMB.annotated.cleaned.tsv.gz")
ok("sQTL duplicate is reported as an sQTL problem",
   grepl("^sQTL fine-mapping table has duplicate", msg_sqtl))
ok("sQTL duplicate is NOT blamed on eQTL or pQTL",
   !grepl("eQTL fine-mapping table", msg_sqtl) && !grepl("pQTL fine-mapping table", msg_sqtl))
ok("a failed coloc_model_input build (sQTL) leaves no cache file behind",
   !file.exists(file.path(DATA_DIR, "coloc_model_input.rds")))

msg_pqtl <- dup_finemap_row("finemapping/coloc/pQTL_susie.COMB.annotated.cleaned.tsv.gz")
ok("pQTL duplicate is reported as a pQTL problem",
   grepl("^pQTL fine-mapping table has duplicate", msg_pqtl))
ok("pQTL duplicate is NOT blamed on eQTL or sQTL",
   !grepl("eQTL fine-mapping table", msg_pqtl) && !grepl("sQTL fine-mapping table", msg_pqtl))
ok("a failed coloc_model_input build (pQTL) leaves no cache file behind",
   !file.exists(file.path(DATA_DIR, "coloc_model_input.rds")))

# Restore the original fixture so anything below (or a re-run) sees the
# well-formed inputs, not the corrupted ones above.
use_fig3_fixture(fx)
source("fig3/00_setup.R")

# --- legitimate data that LOOKS like it could collide must NOT error ------
# Same variant, different genes, and same gene, different variants, are both
# real outcomes (a variant can be fine-mapped for two genes; a gene's
# credible set spans several variants) and must not trip the many-to-one
# guard, which keys on the (variant, gene) PAIR, not either column alone.
fx_legit <- make_fig3_fixture(file.path(tempdir(), "fig3_fixture_legit"))
eqtl_path <- file.path(fx_legit$root, "finemapping/coloc/eQTL_susie.COMB.annotated.cleaned.tsv.gz")
eqtl_legit <- read.delim(eqtl_path)
# Row 1 is (chr1_100_A_G, ENSG001); add the SAME variant for a DIFFERENT gene.
extra <- eqtl_legit[1, ]
extra$molecular_trait_id <- "ENSG002.1"
eqtl_legit <- rbind(eqtl_legit, extra)
write.table(eqtl_legit, eqtl_path, sep = "\t", row.names = FALSE, quote = FALSE)
use_fig3_fixture(fx_legit)
source("fig3/00_setup.R")
mi_legit <- tryCatch(coloc_model_input(refresh = TRUE), error = function(e) e)
ok("same variant, different genes does not error the join",
   !inherits(mi_legit, "error"))

# The base fixture (used for every `mi <- coloc_model_input()` assertion
# above) already has ENSG001 fine-mapped at several distinct eQTL variants --
# that build succeeding is itself proof the same-gene/different-variant case
# does not error; pin it explicitly here too.
ok("same gene, different variants recur in the base fixture",
   length(unique(fm$eQTL$variant[fm$eQTL$molecular_trait_id == "ENSG001"])) > 1)
ok("...and coloc_model_input() built successfully over that fixture",
   is.list(mi) && setequal(names(mi), c("eQTL-sQTL", "eQTL-pQTL", "sQTL-pQTL")))

# Restore the original fixture once more for anything after this.
use_fig3_fixture(fx)
source("fig3/00_setup.R")

cat("test_fig3_builders.R: all assertions passed\n")
