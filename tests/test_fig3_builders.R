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
ok("pQTL keeps every row (`group` is no longer filtered)", nrow(fm$pQTL) == 9L)
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
sq_path <- fixture_path(fx_dup_ann, "susie_sqtl")
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
sig_path <- fixture_path(fx_dup_sig, "enloc_cs", "eQTL-pQTL")
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
# Duplicates on the fine-mapping side are now COLLAPSED by posteriors() rather
# than rejected: intron-level splicing ids make them a routine outcome, not a
# defect (13,950 such pairs in the real sQTL table). The many-to-one guard on the
# posterior join is retained as a cheap invariant -- it now asserts the collapse
# did its job -- but it is no longer reachable from a duplicated input row, so
# the tests that used to drive it have been rewritten around the new contract.
dup_finemap_row <- function(key) {
  fx_dup <- make_fig3_fixture(file.path(tempdir(), paste0("fig3_fixture_dup_", key)))
  p <- fixture_path(fx_dup, key)
  d <- read.delim(p)
  dup <- d[1, ]
  # A DIFFERENT posterior for the same (variant, gene): the collapse must keep
  # the stronger of the two, so a weaker duplicate is not merely tolerated but
  # demonstrably discarded.
  dup$posterior_mean <- -0.01
  write.table(rbind(d, dup), p, sep = "\t", row.names = FALSE, quote = FALSE)
  use_fig3_fixture(fx_dup)
  source("fig3/00_setup.R")
  list(err = tryCatch({ coloc_model_input(refresh = TRUE); NULL }, error = conditionMessage),
       mi  = tryCatch(coloc_model_input(refresh = FALSE), error = function(e) NULL))
}

for (k in c("susie_eqtl", "susie_sqtl", "susie_pqtl")) {
  r <- dup_finemap_row(k)
  ok(paste0(k, ": a duplicate (variant, gene) row no longer fails the build"),
     is.null(r$err))
  if (!is.null(r$mi)) {
    dupd <- vapply(r$mi, function(d)
      if (is.null(d) || !nrow(d)) FALSE else any(duplicated(d[, c("SNP", "Gene")])),
      logical(1))
    ok(paste0(k, ": no (SNP, Gene) pair is duplicated in any pair's model input"),
       !any(dupd))
    weak <- vapply(r$mi, function(d)
      if (is.null(d) || !nrow(d)) FALSE
      else any(abs(d$posterior_mean.x) == 0.01, na.rm = TRUE) ||
           any(abs(d$posterior_mean.y) == 0.01, na.rm = TRUE),
      logical(1))
    ok(paste0(k, ": the weaker duplicate posterior was discarded, not kept"),
       !any(weak))
  }
}

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
eqtl_path <- fixture_path(fx_legit, "susie_eqtl")
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


# =============================================================================
# Shared SuSiE tables with Fig2
# -----------------------------------------------------------------------------
# Fig3 now reads the SAME SuSiE tables as Fig2. Those exports are per-assay COMB
# files, so `group` may be absent, and they are produced by an upstream pipeline
# this repo does not control -- hence an explicit schema check.
# =============================================================================

# Rebuild susie_finemapping() against a fixture mutated by `mutate_fn`, and
# return either the error message or NULL on success.
susie_with <- function(tag, mutate_fn) {
  fx <- make_fig3_fixture(file.path(tempdir(), paste0("fig3_fixture_", tag)))
  mutate_fn(fx)
  use_fig3_fixture(fx)
  source("fig3/00_setup.R")
  tryCatch({ susie_finemapping(refresh = TRUE); NULL }, error = conditionMessage)
}

# `group` is no longer filtered on: these are combined-ancestry exports, one
# file per assay. The real tables hold "COMB.proteomics" etc., so the old
# filter(group == "COMB") would have matched nothing and emptied the table.
fx_grp <- make_fig3_fixture(file.path(tempdir(), "fig3_fixture_grp"))
use_fig3_fixture(fx_grp); source("fig3/00_setup.R")
ok("no row is dropped on `group`",
   nrow(susie_finemapping(refresh = TRUE)$pQTL) ==
     nrow(read.delim(fixture_path(fx_grp, "susie_pqtl"))))

# A missing credible-set column must name every absent column, not just the
# first -- these files come from upstream and a partial message costs a round
# trip per column.
msg2 <- susie_with("nocs", function(fx) {
  p <- fixture_path(fx, "susie_eqtl")
  d <- read.delim(p); d$cs_id <- NULL; d$posterior_sd <- NULL
  write.table(d, p, sep = "\t", row.names = FALSE, quote = FALSE)
})
ok("a SuSiE table missing credible-set columns errors",
   !is.null(msg2))
ok("the error names EVERY missing column, not just the first",
   grepl("cs_id", msg2) && grepl("posterior_sd", msg2))
ok("the error names the offending key and file",
   grepl("susie_eqtl", msg2) && grepl("annotated", msg2))

# =============================================================================
# Duplicate (variant, gene) keys from intron-level splicing ids
# -----------------------------------------------------------------------------
# Splicing molecular_trait_ids are per intron cluster and are collapsed to a
# gene, so one variant can carry a credible set for several introns of the same
# gene (13,950 such pairs in the real sQTL table). The posterior join keeps the
# STRONGEST effect per (variant, gene); keeping both would double-count the
# variant and inflate the regression's N behind panels C and D.
# =============================================================================
fx_dupint <- make_fig3_fixture(file.path(tempdir(), "fig3_fixture_dupint"))
use_fig3_fixture(fx_dupint); source("fig3/00_setup.R")

fm_di <- susie_finemapping(refresh = TRUE)
dup_key <- fm_di$sQTL %>%
  dplyr::count(variant, molecular_trait_id) %>%
  dplyr::filter(n > 1)
ok("the fixture really does carry a duplicate (variant, gene) sQTL key",
   nrow(dup_key) >= 1L)
ok("both intron clusters survive in the fine-mapping table (panels A/B need them)",
   dplyr::n_distinct(fm_di$sQTL$cs_id[fm_di$sQTL$molecular_trait_id == "ENSG001"]) >= 2L)

mi_di <- coloc_model_input(refresh = TRUE)
ok("coloc_model_input() builds despite the duplicate key",
   is.list(mi_di) && length(mi_di) == 3L)

es <- mi_di[["eQTL-sQTL"]]
ok("no (SNP, Gene) pair is duplicated after the posterior join",
   !any(duplicated(es[, c("SNP", "Gene")])))
# Only the DUPLICATED (variant, gene) pair should show the stronger posterior;
# ENSG001's other variants keep their own, much smaller, effects.
row <- es[es$SNP == "chr1_100_A_G" & es$Gene == "ENSG001", ]
ok("the duplicated (variant, gene) pair appears exactly once", nrow(row) == 1L)
ok("the STRONGEST of its two intron clusters is the one retained",
   nrow(row) == 1L && isTRUE(abs(abs(row$posterior_mean.y[1]) - 0.95) < 1e-8))
cat("test_fig3_builders.R: all assertions passed\n")
