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
ok("carries the genuinely variant-level annotation columns",
   all(c("consequence", "phylop") %in% names(map)))
ok("does NOT carry the gene-relative columns -- those are variant x gene, not variant-level",
   !any(c("distTSS", "gene_type") %in% names(map)))
ok("carries a per-variant median credible-set size", "cs_size_med" %in% names(map))
ok("cs_size_med is numeric and positive",
   is.numeric(map$cs_size_med) && all(map$cs_size_med > 0, na.rm = TRUE))
ok("allele frequencies survive the widening",
   all(AF_POP_COLS %in% names(map)))

# --- gene-relative feature map: keyed on (variant, gene), not variant alone -
# The fixture adds a second fine-mapped occurrence of variant 5
# (snps[5] = "chr1_500_A_G") under ENSG002, with a DIFFERENT gene-relative
# annotation (distTSS = 50000, gene_type = "lncRNA") than it carries for
# ENSG001 (distTSS = 500, gene_type = "protein_coding") -- see
# tests/fixtures.R. Building .variant_gene_feature_map() must not error (the
# defect this task fixes: the ORIGINAL all-columns guard hard-errored on
# exactly this shape) and must keep both rows, one per gene.
gmap <- .variant_gene_feature_map(fm)

ok(".variant_gene_feature_map() does not error on a variant shared by two genes with differing gene-relative values",
   is.data.frame(gmap))
ok("one row per (variant, gene) pair, not per variant",
   !any(duplicated(gmap[c("variant", "gene")])))

dual_rows <- gmap[gmap$variant == "chr1_500_A_G", ]
ok("the shared variant appears under BOTH genes",
   setequal(dual_rows$gene, c("ENSG001", "ENSG002")))
ok("distTSS differs PER GENE for the shared variant",
   dual_rows$distTSS[dual_rows$gene == "ENSG001"] == 500 &&
   dual_rows$distTSS[dual_rows$gene == "ENSG002"] == 50000)
ok("gene_type differs PER GENE for the shared variant",
   dual_rows$gene_type[dual_rows$gene == "ENSG001"] == "protein_coding" &&
   dual_rows$gene_type[dual_rows$gene == "ENSG002"] == "lncRNA")

# --- tss_proximal / protein_coding are computed PER GENE via coloc_variant_level() ---
# This is the assertion that proves the fix is real rather than just a
# suppressed guard: the SAME variant, joined against two different Gene
# rows in the enloc variant-level table, must come out TSS-proximal +
# protein-coding for one gene and neither for the other.
vlevel <- coloc_variant_level(refresh = TRUE)
# Restricted to the eQTL-pQTL pair, the only one carrying an ENSG002 row for
# this variant (see the dual_gene_row appended to that pair's enloc_variant
# fixture above) -- the other two pairs only ever see it under ENSG001, so
# without this filter there would be three ENSG001 rows to one ENSG002 row.
dual_vlevel <- vlevel[vlevel$SNP == "chr1_500_A_G" & vlevel$type == "eQTL-pQTL" &
                        vlevel$Gene %in% c("ENSG001", "ENSG002"), ]

ok("coloc_variant_level() carries both (variant, gene) rows for the shared variant",
   setequal(dual_vlevel$Gene, c("ENSG001", "ENSG002")))

# build_feature_matrix() needs max_af, which coloc_feature_enrichment() folds
# in before calling it (see Fig3/00_setup.R); replicate that fold here rather
# than pull the memoized coloc_feature_enrichment() builder, so this stays a
# direct, minimal exercise of build_feature_matrix() on the two dual-gene rows.
dual_folded <- dual_vlevel[AF_POP_COLS]
dual_folded[] <- lapply(dual_folded, function(x) ifelse(x > .5, 1 - x, x))
dual_vlevel$max_af <- do.call(pmax, as.list(dual_folded))

dual_feat <- build_feature_matrix(dual_vlevel)
ok("tss_proximal is TRUE for ENSG001 and FALSE for ENSG002, for the SAME variant",
   isTRUE(dual_feat$tss_proximal[dual_feat$Gene == "ENSG001"]) &&
   isFALSE(dual_feat$tss_proximal[dual_feat$Gene == "ENSG002"]))
ok("protein_coding is TRUE for ENSG001 and FALSE for ENSG002, for the SAME variant",
   isTRUE(dual_feat$protein_coding[dual_feat$Gene == "ENSG001"]) &&
   isFALSE(dual_feat$protein_coding[dual_feat$Gene == "ENSG002"]))

# --- the variant-level guard still fires on a genuinely variant-level conflict ---
# A disagreement in an AF column (variant-level, not gene-relative) across
# assay tables must still hard-error -- this is the guard
# tests/test_fig3_builders.R already exercises for coloc_variant_level(); this
# repeats it directly against .variant_feature_map() so it stays pinned next
# to the gene-relative split above.
fx_dup_af <- make_fig3_fixture(file.path(tempdir(), "fig3_enrich_dup_af"))
sq_dup_path <- fixture_path(fx_dup_af, "susie_sqtl")
sq_dup <- read.delim(sq_dup_path)
sq_dup$gvs_afr_af[1] <- sq_dup$gvs_afr_af[1] + 0.321   # disagree with eQTL/pQTL
write.table(sq_dup, sq_dup_path, sep = "\t", row.names = FALSE, quote = FALSE)
use_fig3_fixture(fx_dup_af)
source("fig3/00_setup.R")
expect_error("a genuinely variant-level column (AF) disagreeing across assays still errors .variant_feature_map()",
             .variant_feature_map(susie_finemapping(refresh = TRUE)))

# --- intra-table (variant, gene) duplicates are collapsed, cross-table ones still error ---
# Real-data finding (see NOTES.md's "Duplicate (variant, gene) keys are
# collapsed, not rejected"): a SINGLE SuSiE table can itself carry two rows
# for the same (variant, gene) -- an upstream ENCODE dual-cCRE-match
# artifact -- that disagree on a variant-level cCRE flag. Collapsing those
# (first row wins) must not also swallow a genuine disagreement BETWEEN two
# DIFFERENT assay tables about the same (variant, gene), which is exactly
# what the guard exists to catch. In-memory frames built by hand (not via
# the file-based fixture) so this stays a minimal, direct exercise of
# .variant_feature_map() itself.
mk_af_row <- function(variant, molecular_trait_id, ctcf_bound, af = 0.05) {
  af_row    <- stats::setNames(as.list(rep(af, length(AF_POP_COLS))), AF_POP_COLS)
  ccre_row  <- stats::setNames(as.list(rep(FALSE, length(FEATURE_CCRE_COLS))), FEATURE_CCRE_COLS)
  ccre_row[["CTCF.bound"]] <- ctcf_bound
  tibble::as_tibble(c(
    list(variant = variant, molecular_trait_id = molecular_trait_id),
    af_row, ccre_row,
    list(consequence = "intron_variant", phylop = 1.0, cs_size = 3L)))
}
empty_af_table <- mk_af_row("chrZ_0_A_G", "ENSGZ", ctcf_bound = FALSE)[0, ]

intratab <- dplyr::bind_rows(
  mk_af_row("chrX_1_A_G", "ENSGX", ctcf_bound = TRUE),
  mk_af_row("chrX_1_A_G", "ENSGX", ctcf_bound = FALSE))    # same table, same (variant, gene)
fm_intratab <- list(eQTL = intratab, sQTL = empty_af_table, pQTL = empty_af_table)
map_intratab <- tryCatch(.variant_feature_map(fm_intratab), error = function(e) e)
ok("an intra-table duplicate (variant, gene) row is collapsed (first row wins), not rejected",
   !inherits(map_intratab, "error") && nrow(map_intratab) == 1L)

fm_crosstab <- list(
  eQTL = mk_af_row("chrX_2_A_G", "ENSGX", ctcf_bound = TRUE),
  sQTL = mk_af_row("chrX_2_A_G", "ENSGX", ctcf_bound = FALSE),  # DIFFERENT table, same (variant, gene)
  pQTL = empty_af_table)
err_crosstab <- tryCatch({ .variant_feature_map(fm_crosstab); NULL }, error = conditionMessage)
ok("the SAME disagreement BETWEEN two different assay tables still errors -- the collapse is scoped within one table only",
   !is.null(err_crosstab))

# Restore the original enrichment fixture for anything appended below.
use_fig3_fixture(fx)
source("fig3/00_setup.R")

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

# --- per-feature model fit --------------------------------------------------
set.seed(2026)
n <- 800
sim <- data.frame(
  PIP_qtl     = runif(n, 0.1, 0.9),
  cs_size_med = sample(2:40, n, replace = TRUE),
  # `hit` is strongly associated with the outcome; `dud` is not.
  hit         = rep(c(TRUE, FALSE), each = n / 2),
  dud         = sample(c(TRUE, FALSE), n, replace = TRUE),
  # `thin` has only 4 carriers -- below the (current default) floor of 5.
  thin        = c(rep(TRUE, 4), rep(FALSE, n - 4)))
sim$outcome <- c(rbinom(n / 2, 1, 0.65), rbinom(n / 2, 1, 0.20)) == 1L

res <- fit_feature_enrichment(sim, "outcome", c("hit", "dud", "thin"))

ok("returns one row per surviving feature", nrow(res) == 2L)
ok("a below-floor feature is dropped", !"thin" %in% res$feature)
ok("an enriched feature has OR > 1 with a CI excluding 1",
   res$odds_ratio[res$feature == "hit"] > 1 && res$conf.low[res$feature == "hit"] > 1)
ok("a null feature's CI includes 1",
   res$conf.low[res$feature == "dud"] < 1 & res$conf.high[res$feature == "dud"] > 1)
ok("n_carriers is reported", all(res$n_carriers > 0))
ok("the tidy shape is exactly as documented",
   identical(sort(names(res)),
             sort(c("feature", "odds_ratio", "conf.low", "conf.high",
                    "p.value", "n_carriers"))))

msg <- capture.output(
  fit_feature_enrichment(sim, "outcome", c("hit", "thin")), type = "message")
ok("the dropped feature is NAMED in a message, not silently discarded",
   any(grepl("thin", msg)))

ok("no surviving feature returns a zero-row frame of the right shape",
   nrow(fit_feature_enrichment(sim, "outcome", "thin")) == 0L)

# --- a per-feature fit/profile failure is isolated, not fatal --------------
# `separated` clears the carrier floor (its smallest table(x, outcome) cell
# is well above the current default of 5) but is constructed so that
# INCLUDING it as a glm
# covariate makes the outcome perfectly linearly separable: rows with
# `separated == TRUE` flip to TRUE past a PIP_qtl of 0.3, rows with
# `separated == FALSE` flip past 0.7. An additive glm term represents that
# as an intercept shift, so `outcome ~ separated + PIP_qtl + log10(cs_size)`
# is a clean separating hyperplane and stats::confint() dies with "need at
# least two non-NA values to interpolate" -- the reviewer's reproduction of
# the real failure mode. `healthy` has no such construction: with
# `separated` excluded from ITS model, the same PIP_qtl value can carry
# either outcome depending on the latent (unmodelled) `separated` status, so
# `healthy`'s own fit is well-behaved and must still be returned.
set.seed(2026)
n3 <- 100
sep_frame <- data.frame(
  PIP_qtl     = runif(n3, 0, 1),
  cs_size_med = sample(2:40, n3, replace = TRUE),
  separated   = c(rep(TRUE, 40L), rep(FALSE, 60L)),
  healthy     = sample(c(TRUE, FALSE), n3, replace = TRUE))
sep_frame$outcome <- with(sep_frame, ifelse(separated, PIP_qtl > 0.3, PIP_qtl > 0.7))

msg2 <- capture.output(
  res2 <- fit_feature_enrichment(sep_frame, "outcome", c("separated", "healthy")),
  type = "message")

ok("a fit/profile failure on one feature does not error the whole call",
   is.data.frame(res2))
ok("the healthy feature is still returned despite the other's failure",
   "healthy" %in% res2$feature)
ok("the separated (unfittable) feature is absent from the results",
   !"separated" %in% res2$feature)
ok("the separated feature's failure is NAMED in a message, not swallowed",
   any(grepl("separated", msg2)))

# --- builder ----------------------------------------------------------------
enr <- coloc_feature_enrichment(refresh = TRUE)

# The 9-row Fig3 fixture (see tests/fixtures.R) is sized for panels C/D's
# per-MAF-bin regressions, not for meaningfully covering the carrier floor
# within a single assay pair. At the OLD floor of 10, min(table(feature,
# outcome)) could never reach 10 out of 9 total rows, so nothing ever
# survived. At the current default of 5 it clears for exactly one
# (analysis, pair, feature) combination -- analysis B / eQTL-pQTL /
# tss_proximal, n_carriers = 5 -- and nothing for analysis A, so `enr` here
# is a well-formed 1-row frame, not the fully-empty one the old floor
# produced. The assertions below therefore check the SHAPE of the result
# unconditionally, and check per-row properties only if any rows happen to
# exist -- exactly per the "do not depend on every combination succeeding"
# contract; they do NOT assume both analyses necessarily appear on a
# fixture this small (the `.enrich_pairs()` synthetic-matrix block below
# covers that property on a properly-sized input instead).
ok("the documented columns are all present",
   all(c("analysis", "pair", "feature", "odds_ratio", "conf.low", "conf.high",
         "p.value", "padj", "n_carriers", "flagged") %in% names(enr)))
ok("coloc_feature_enrichment() returns a well-formed tibble even when the fixture is too small for full coverage",
   is.data.frame(enr))
ok("only the documented analysis labels ever appear",
   all(unique(enr$analysis) %in% c("A", "B")))
ok("both analyses return the SAME columns",
   identical(sort(names(enr[enr$analysis == "A", ])),
             sort(names(enr[enr$analysis == "B", ]))))

if (nrow(enr)) {
  ok("padj is FDR within (analysis, pair), so never below the raw p",
     all(enr$padj >= enr$p.value - 1e-9))
  ok("flagged is logical", is.logical(enr$flagged))
  ok("odds ratios are finite -- the carrier floor rules out separation",
     all(is.finite(enr$odds_ratio)))
} else {
  message("  note: the fixture is too small for any (analysis, pair) ",
          "combination to clear the 5-carrier floor; coloc_feature_enrichment() ",
          "still returned a well-formed zero-row tibble of the documented shape.")
  ok("flagged column is logical-typed even with zero rows", is.logical(enr$flagged))
}

# --- .enrich_pairs() orchestration on a synthetic feature matrix ------------
# The real Fig3 fixture (9 rows) barely clears the current default 5-carrier
# floor (exactly 1 row, analysis B only -- see the builder block above), so
# it is still far too thin to exercise the orchestration for real: deleting
# the FDR line, reverting `flagged`'s denominator, or dropping the `rare`
# exclusion from B could all still leave `enr` looking plausible on a single
# row. This block builds a small, fully
# deterministic feature matrix by hand (in the spirit of `sim` above) and
# calls the pure `.enrich_pairs(v, min_carriers = ...)` helper directly --
# bypassing coloc_variant_level()/build_feature_matrix() and the memoized
# cache entirely -- so the carrier floor can be lowered without touching
# tests/fixtures.R and every row asserted on is a REAL, fitted result.
#
# Two pairs, each split 50/50 on Coloc, n_half rows per Coloc value:
#   PAIR_OK  -- `hit` is a real, strongly-associated feature (minority cell
#               size 15, comfortably above the floor); `dud1`/`dud2`/`dud3`
#               and `rare` are balanced-but-unassociated (~30/30 cells) --
#               ALL 5 candidates clear the floor in analysis A.
#   PAIR_LOW -- same `hit` column (so its fit is identical to PAIR_OK's,
#               which is what makes the pooled-vs-within-group FDR
#               divergence below exact and reproducible), but `dud1`/`dud2`/
#               `dud3`/`rare` are made CONSTANT columns (zero variance
#               against Coloc), so only `hit` clears the floor: 1 of 5
#               candidates survives in analysis A, well under half.
mk_enrich_pair <- function(pair, n_half, minority, sparse) {
  Coloc <- c(rep(TRUE, n_half), rep(FALSE, n_half))
  # `hit`: minority-cell size is fixed at `minority` in BOTH Coloc groups, so
  # table(hit, Coloc) has all four cells >= minority regardless of `sparse`.
  hit <- c(rep(TRUE, n_half - minority), rep(FALSE, minority),
           rep(FALSE, n_half - minority), rep(TRUE, minority))
  if (!sparse) {
    # Period-2/4 alternating patterns: balanced against Coloc (~n_half/2 per
    # cell) and, restricted to the Coloc == TRUE half used by analysis B,
    # balanced against `rare` too -- every candidate clears the floor in
    # both analyses.
    dud1 <- rep(c(TRUE, FALSE),             n_half)
    dud2 <- rep(c(FALSE, TRUE),             n_half)
    dud3 <- rep(c(TRUE, TRUE, FALSE, FALSE), n_half / 2)
    rare <- rep(c(TRUE, FALSE, FALSE, TRUE), n_half / 2)
  } else {
    # Zero variance: table(x, Coloc) has an all-zero row/column no matter how
    # the floor is set, so these four are dropped unconditionally.
    dud1 <- rep(FALSE, 2 * n_half)
    dud2 <- rep(FALSE, 2 * n_half)
    dud3 <- rep(TRUE,  2 * n_half)
    rare <- rep(FALSE, 2 * n_half)
  }
  # Analysis B now regresses log10(max_af) directly rather than the binary
  # `rare` outcome, so the synthetic matrix needs a continuous, strictly
  # positive max_af. A smooth deterministic ramp (independent of `rare`,
  # which keeps its own hand-crafted pattern above for analysis A's floor/
  # sparse behaviour) gives every feature's fit real variation to work with
  # without introducing randomness into an otherwise fully deterministic
  # fixture.
  max_af <- rep(seq(0.01, 0.35, length.out = n_half), 2)
  tibble::tibble(
    type = pair, Coloc = Coloc, hit = hit, dud1 = dud1, dud2 = dud2, dud3 = dud3,
    rare = rare, max_af = max_af,
    PIP_qtl     = rep(seq(0.1, 0.9, length.out = n_half), 2),
    cs_size_med = rep(2:(n_half + 1), 2)[seq_len(2 * n_half)])
}

v_synth <- dplyr::bind_rows(
  mk_enrich_pair("PAIR_OK",  n_half = 60, minority = 15, sparse = FALSE),
  mk_enrich_pair("PAIR_LOW", n_half = 60, minority = 15, sparse = TRUE))

enr2 <- .enrich_pairs(v_synth, min_carriers = 5)

ok("candidate features clear the floor on the synthetic matrix (not vacuous: rows exist)",
   nrow(enr2) > 0L)
ok("both analyses are present on the synthetic matrix",
   setequal(unique(enr2$analysis), c("A", "B")))
ok("rare is a tested feature in analysis A",
   "rare" %in% enr2$feature[enr2$analysis == "A"])
ok("rare NEVER appears among analysis B's features -- it is B's outcome",
   !"rare" %in% enr2$feature[enr2$analysis == "B"])

ok("the underpowered pair (1/5 candidates survive) is flagged",
   isTRUE(enr2$flagged[enr2$analysis == "A" & enr2$pair == "PAIR_LOW" & enr2$feature == "hit"]))
ok("the healthy pair (5/5 candidates survive) is NOT flagged",
   all(!enr2$flagged[enr2$analysis == "A" & enr2$pair == "PAIR_OK"]))

# FDR-within-(analysis, pair) vs. a pooled correction, checked on a
# MULTI-feature group (analysis B / PAIR_OK has 4 surviving candidates: hit,
# dud1, dud2, dud3). A single-feature group is the WRONG group to check this
# on: BH applied to one p-value is always a no-op (padj == p.value) whether
# or not p.adjust() actually runs, so such a check would pass even with the
# p.adjust() line deleted. A 4-element group does not have that blind spot.
grp_b_ok <- which(enr2$analysis == "B" & enr2$pair == "PAIR_OK")
stopifnot(length(grp_b_ok) == 4L)
expected_within_b_ok <- stats::p.adjust(enr2$p.value[grp_b_ok], method = "fdr")
pooled_padj_all       <- stats::p.adjust(enr2$p.value, method = "fdr")

ok("padj on PAIR_OK's analysis-B group matches the WITHIN-(analysis,pair) BH adjustment, feature by feature",
   isTRUE(all.equal(enr2$padj[grp_b_ok], expected_within_b_ok)))
ok("at least one of those padj values differs from its own raw p-value -- a real correction is applied, not a p.value passthrough",
   any(abs(enr2$padj[grp_b_ok] - enr2$p.value[grp_b_ok]) > 1e-9))
ok("at least one of those padj values differs from a pooled-across-everything BH correction -- proving padj is scoped per (analysis, pair), not pooled",
   any(abs(enr2$padj[grp_b_ok] - pooled_padj_all[grp_b_ok]) > 1e-9))

ok("every returned odds ratio on the synthetic matrix is finite",
   all(is.finite(enr2$odds_ratio)))

# --- panel E -----------------------------------------------------------------
pE <- build_panelE(enr)
ok("panel E is a ggplot", inherits(pE, "ggplot"))
ok("panel E plots analysis A only", all(pE$data$analysis == "A"))
ok("panel E is a two-pair contrast, one row per feature",
   !any(duplicated(pE$data$feature)))

# The 9-row Fig3 fixture makes `enr` (built above) a well-formed but
# ZERO-ROW tibble, so the three assertions just above are trivially true of
# an empty data frame and cannot catch a broken pivot, a wrong analysis
# filter, or a missing-pair leak. A small synthetic enrichment frame -- same
# schema as `enr`/`enr2` above, built by hand rather than through the
# memoized pipeline -- exercises build_panelE() on REAL rows instead.
enr_panelE <- tibble::tibble(
  analysis = c("A", "A",             # f1: both required pairs -> kept
               "A",                  # f2: eQTL-pQTL only       -> dropped
               "A",                  # f3: eQTL-sQTL only       -> dropped
               "A", "A", "A",        # f4: both pairs + the structurally-
                                      #     absent third pair    -> kept,
                                      #     third pair ignored
               "B", "B"),            # f5: analysis B           -> dropped
  pair     = c("eQTL-pQTL", "eQTL-sQTL",
               "eQTL-pQTL",
               "eQTL-sQTL",
               "eQTL-pQTL", "eQTL-sQTL", "sQTL-pQTL",
               "eQTL-pQTL", "eQTL-sQTL"),
  feature    = c("f1", "f1", "f2", "f3", "f4", "f4", "f4", "f5", "f5"),
  odds_ratio = c(2.0, 1.5, 3.0, 0.5, 1.2, 0.8, 4.0, 2.0, 2.0),
  padj       = c(0.01, 0.20, 0.01, 0.01, 0.04, 0.04, 0.01, 0.01, 0.01))

pE2 <- build_panelE(enr_panelE)

ok("synthetic panel E: analysis A only, B is excluded",
   all(pE2$data$analysis == "A"))
ok("synthetic panel E: one row per feature",
   !any(duplicated(pE2$data$feature)))
ok("synthetic panel E: both required pairs are pivoted into columns",
   all(c("log_or_eQTL-pQTL", "log_or_eQTL-sQTL") %in% names(pE2$data)))
ok("synthetic panel E: the structurally-absent sQTL-pQTL pair is not pivoted in",
   !any(grepl("sQTL-pQTL", names(pE2$data), fixed = TRUE)))
ok("synthetic panel E: features missing from either pair are excluded",
   setequal(pE2$data$feature, c("f1", "f4")))
ok("synthetic panel E: log-OR is log2 of the odds ratio",
   isTRUE(all.equal(pE2$data$`log_or_eQTL-pQTL`[pE2$data$feature == "f1"],
                    log2(2.0))))

# --- enrichment_forest_panel() on the zero-row adversarial cases -----------
# The supplement's forest() (now enrichment_forest_panel() in Fig3/00_setup.R,
# extracted specifically so this is testable) used to facet_wrap(~pair)
# unconditionally, which errors ("Faceting variables must have at least one
# value") at BUILD/render time -- not at ggplot() construction -- whenever
# its input has zero rows. That happens not only when coloc_feature_enrichment()
# is fully empty (the 9-row fixture's real, documented outcome above) but also
# when only ONE analysis half is empty while the other has real rows -- a
# state .enrich_pairs() explicitly allows. ggplot2::ggplot_build() is called
# below to force evaluation of the facet spec, matching what ggsave() would
# trigger; a bare ggplot() call is lazy and would NOT reproduce the crash.
schema_cols <- c("analysis", "pair", "feature", "odds_ratio", "conf.low",
                 "conf.high", "p.value", "padj", "n_carriers", "flagged")

empty_enr <- tibble::tibble(
  analysis = character(), pair = character(), feature = character(),
  odds_ratio = numeric(), conf.low = numeric(), conf.high = numeric(),
  p.value = numeric(), padj = numeric(), n_carriers = integer(),
  flagged = logical())

mixed_enr <- tibble::tibble(
  analysis   = c("A", "A"),
  pair       = c("eQTL-pQTL", "eQTL-sQTL"),
  feature    = c("f1", "f2"),
  odds_ratio = c(2.0, 1.5),
  conf.low   = c(1.1, 0.9),
  conf.high  = c(3.5, 2.4),
  p.value    = c(0.01, 0.2),
  padj       = c(0.01, 0.2),
  n_carriers = c(20L, 25L),
  flagged    = c(FALSE, FALSE))

ok("the synthetic empty frame carries the documented schema",
   identical(sort(names(empty_enr)), sort(schema_cols)))
ok("the synthetic mixed frame carries the documented schema",
   identical(sort(names(mixed_enr)), sort(schema_cols)))

# Case (a): fully empty input (reviewer's first reproduction).
p_empty <- try(
  ggplot2::ggplot_build(
    enrichment_forest_panel(dplyr::filter(empty_enr, analysis == "A"), "A")),
  silent = TRUE)
ok("a fully empty analysis frame does not error on build (was: 'Faceting variables must have at least one value')",
   !inherits(p_empty, "try-error"))
ok("enrichment_forest_panel() on a fully empty frame still returns a ggplot object",
   inherits(enrichment_forest_panel(dplyr::filter(empty_enr, analysis == "A"), "A"), "ggplot"))

# Case (b): analysis A has real rows, analysis B is empty (reviewer's second,
# adversarial reproduction -- NOT fully empty overall).
p_b_empty <- try(
  ggplot2::ggplot_build(
    enrichment_forest_panel(dplyr::filter(mixed_enr, analysis == "B"), "B")),
  silent = TRUE)
ok("analysis B empty while analysis A is non-empty does not error on build",
   !inherits(p_b_empty, "try-error"))
p_b_panel <- enrichment_forest_panel(dplyr::filter(mixed_enr, analysis == "B"), "B")
ok("the empty-B placeholder is still a ggplot object",
   inherits(p_b_panel, "ggplot"))
ok("the empty-B placeholder names the analysis so it reads as 'no results', not 'not run'",
   any(grepl("B", vapply(p_b_panel$layers, function(l)
     if (!is.null(l$aes_params$label)) l$aes_params$label else NA_character_,
     character(1)), fixed = TRUE)))

# The non-empty half (analysis A here) must still build normally -- the guard
# must not degrade the populated case into a placeholder too.
p_a_populated <- try(
  ggplot2::ggplot_build(
    enrichment_forest_panel(dplyr::filter(mixed_enr, analysis == "A"), "A")),
  silent = TRUE)
ok("the non-empty half still builds as a real faceted forest plot, not a placeholder",
   !inherits(p_a_populated, "try-error"))

# --- analysis B is continuous log MAF ---------------------------------------
# Build a matrix where `hi_maf` carriers are systematically COMMONER and
# `lo_maf` carriers systematically RARER among colocalizing variants.
set.seed(4242)
nB <- 400
mB <- data.frame(
  type        = "eQTL-pQTL",
  Coloc       = TRUE,
  PIP_qtl     = runif(nB, .1, .9),
  cs_size_med = sample(2:30, nB, replace = TRUE),
  hi_maf      = rep(c(TRUE, FALSE), each = nB / 2),
  lo_maf      = rep(c(FALSE, TRUE), each = nB / 2))
# carriers of hi_maf sit near MAF 30%, non-carriers near 3%
mB$max_af <- ifelse(mB$hi_maf, runif(nB, .2, .45), runif(nB, .01, .05))
mB$rare   <- mB$max_af < .01

resB <- .enrich_pairs(mB, min_carriers = 10)
b <- resB[resB$analysis == "B", ]

ok("analysis B returns rows on a continuous outcome", nrow(b) > 0L)
ok("analysis B is labelled as a MAF fold change",
   all(b$estimate_type[b$analysis == "B"] == "maf_fold_change"))
ok("analysis A is still an odds ratio",
   all(resB$estimate_type[resB$analysis == "A"] == "odds_ratio"))
ok("a feature whose carriers are COMMONER gives fold-change < 1",
   b$odds_ratio[b$feature == "hi_maf"] < 1)
ok("a feature whose carriers are RARER gives fold-change > 1",
   b$odds_ratio[b$feature == "lo_maf"] > 1)
ok("the CI brackets the estimate on the same scale",
   all(b$conf.low < b$odds_ratio & b$odds_ratio < b$conf.high, na.rm = TRUE))

# --- coefficient-row lookup handles digit-leading (non-syntactic) names ----
# R keeps backticks around a formula term whose name is not syntactically
# valid (e.g. anything starting with a digit), so its glm coefficient row is
# named `` `3_prime_UTR_variant`TRUE ``, not the naive "3_prime_UTR_variantTRUE".
# The lookup in fit_feature_enrichment() must resolve both forms, exactly --
# not via a substring/grep() match, since a feature name can be a substring of
# another model term ("qtl" inside "PIP_qtl"). This block pins that directly:
# a digit-leading feature must be FITTED AND RETURNED (not silently dropped as
# a "failed to fit/profile" case), alongside an ordinary feature in the same
# call, under both families used by the two analyses (binomial for analysis A,
# gaussian for analysis B). The digit-leading feature is given a strong, known
# effect and the ordinary feature a null effect so a lookup that silently grabs
# the WRONG row (estimates swapped) is caught, not just a dropped-feature bug.
set.seed(99)
n4 <- 800
digit_frame <- data.frame(
  PIP_qtl     = runif(n4, 0.1, 0.9),
  cs_size_med = sample(2:40, n4, replace = TRUE),
  # Strong, known positive association with the binomial outcome.
  `3_prime_UTR_variant` = rep(c(TRUE, FALSE), each = n4 / 2),
  # Null / no association -- an ordinary, syntactically valid name.
  PLS         = sample(c(TRUE, FALSE), n4, replace = TRUE),
  check.names = FALSE)
digit_frame$outcome <- c(rbinom(n4 / 2, 1, 0.75), rbinom(n4 / 2, 1, 0.20)) == 1L

msg_digit <- capture.output(
  res_digit <- fit_feature_enrichment(digit_frame, "outcome",
                                      c("3_prime_UTR_variant", "PLS")),
  type = "message")

ok("a digit-leading feature name is fitted and returned, not dropped (binomial)",
   "3_prime_UTR_variant" %in% res_digit$feature)
ok("no 'failed to fit/profile' message names the digit-leading feature (binomial)",
   !any(grepl("3_prime_UTR_variant", msg_digit) & grepl("failed to fit", msg_digit)))
ok("an ordinary feature name still works alongside it in the same call (binomial)",
   "PLS" %in% res_digit$feature)
ok("the digit-leading feature's estimate reflects its strong effect, not swapped with PLS's (binomial)",
   res_digit$odds_ratio[res_digit$feature == "3_prime_UTR_variant"] > 1 &&
   res_digit$conf.low[res_digit$feature == "3_prime_UTR_variant"] > 1)
ok("the ordinary (null) feature's CI still includes 1, unaffected by the lookup fix (binomial)",
   res_digit$conf.low[res_digit$feature == "PLS"] < 1 &&
   res_digit$conf.high[res_digit$feature == "PLS"] > 1)

# Same shape under family = "gaussian" (analysis B's path): outcome is
# continuous, digit-leading feature carries a strong, known effect on it.
set.seed(100)
digit_frame_g <- data.frame(
  PIP_qtl     = runif(n4, 0.1, 0.9),
  cs_size_med = sample(2:40, n4, replace = TRUE),
  `3_prime_UTR_variant` = rep(c(TRUE, FALSE), each = n4 / 2),
  PLS         = sample(c(TRUE, FALSE), n4, replace = TRUE),
  check.names = FALSE)
digit_frame_g$outcome <- ifelse(digit_frame_g$`3_prime_UTR_variant`,
                                rnorm(n4, 5, 1), rnorm(n4, 0, 1))

msg_digit_g <- capture.output(
  res_digit_g <- fit_feature_enrichment(digit_frame_g, "outcome",
                                        c("3_prime_UTR_variant", "PLS"),
                                        family = "gaussian"),
  type = "message")

ok("a digit-leading feature name is fitted and returned, not dropped (gaussian)",
   "3_prime_UTR_variant" %in% res_digit_g$feature)
ok("no 'failed to fit/profile' message names the digit-leading feature (gaussian)",
   !any(grepl("3_prime_UTR_variant", msg_digit_g) & grepl("failed to fit", msg_digit_g)))
ok("an ordinary feature name still works alongside it in the same call (gaussian)",
   "PLS" %in% res_digit_g$feature)
ok("the digit-leading feature's coefficient reflects its known strong effect, not swapped with PLS's (gaussian)",
   abs(res_digit_g$odds_ratio[res_digit_g$feature == "3_prime_UTR_variant"] -
       res_digit_g$odds_ratio[res_digit_g$feature == "PLS"]) > 0.5)

# --- the carrier floor is 5, not 10 -----------------------------------------
# Measured on the real staged data (see NOTES.md and fit_feature_enrichment()'s
# comment): floor 10 -> 34 FDR-significant features, floor 5 -> 40, and every
# floor-10 significant feature survives at floor 5. These tests pin the two
# edges of that change directly against fit_feature_enrichment(), independent
# of the real data.

# A feature with 5-9 total carriers: dropped at the OLD default (10), tested
# at the CURRENT default (5). Uses family = "gaussian" (analysis B's path) so
# the floor applies directly to table(feature) -- i.e. the total carrier
# count -- rather than to a 2x2 cell that a total-carrier count alone does not
# pin down.
set.seed(3030)
n5 <- 200
mid_frame <- data.frame(
  PIP_qtl     = runif(n5, 0.1, 0.9),
  cs_size_med = sample(2:40, n5, replace = TRUE),
  # `mid` has exactly 7 carriers: inside [5, 9], the band this task adds
  # coverage for.
  mid         = c(rep(TRUE, 7), rep(FALSE, n5 - 7)))
mid_frame$outcome <- rnorm(n5, ifelse(mid_frame$mid, 5, 0), 1)

res_mid_default <- fit_feature_enrichment(mid_frame, "outcome", "mid", family = "gaussian")
ok("a feature with 7 carriers (5-9 band) is TESTED at the current default floor of 5",
   "mid" %in% res_mid_default$feature)
ok("its reported n_carriers is 7",
   res_mid_default$n_carriers[res_mid_default$feature == "mid"] == 7L)

msg_mid_old <- capture.output(
  res_mid_old <- fit_feature_enrichment(mid_frame, "outcome", "mid",
                                        family = "gaussian", min_carriers = 10),
  type = "message")
ok("the SAME feature would have been dropped under the old default of 10",
   !"mid" %in% res_mid_old$feature)
ok("...and named in the dropped-feature message under the old floor",
   any(grepl("mid", msg_mid_old)))

# A feature with fewer than 5 carriers is still dropped at the current
# default, and still named in the message -- the floor moved, it did not
# disappear.
set.seed(3031)
verylow_frame <- data.frame(
  PIP_qtl     = runif(n5, 0.1, 0.9),
  cs_size_med = sample(2:40, n5, replace = TRUE),
  verylow     = c(rep(TRUE, 3), rep(FALSE, n5 - 3)))
verylow_frame$outcome <- rnorm(n5, ifelse(verylow_frame$verylow, 5, 0), 1)

msg_verylow <- capture.output(
  res_verylow <- fit_feature_enrichment(verylow_frame, "outcome", "verylow",
                                        family = "gaussian"),
  type = "message")
ok("a feature with 3 carriers (below the current floor of 5) is still dropped",
   !"verylow" %in% res_verylow$feature)
ok("...and is still named in the dropped-feature message",
   any(grepl("verylow", msg_verylow)))

# --- low-n points render distinctly but do not disappear -------------------
# build_panelE() and enrichment_forest_panel() must both keep rendering when
# their input mixes low-n (n_carriers < 10, a lead) and normal-n rows -- the
# whole point of task 2 is that low-n results stay ON the plot, visually
# weaker, rather than being filtered out like the carrier floor filters out
# sub-floor features upstream.
mix_panelE <- tibble::tibble(
  analysis   = rep("A", 4),
  pair       = rep(c("eQTL-pQTL", "eQTL-sQTL"), 2),
  feature    = rep(c("f_lo", "f_hi"), each = 2),
  odds_ratio = c(2.0, 1.8, 3.0, 2.5),
  padj       = c(0.01, 0.01, 0.02, 0.02),
  # f_lo: min(6, 20) = 6 < 10 -> a lead. f_hi: min(25, 30) = 25 >= 10 -> normal.
  n_carriers = c(6L, 20L, 25L, 30L))

p_mixE <- build_panelE(mix_panelE)
ok("build_panelE() returns a ggplot object on a mix of low-n and normal-n rows",
   inherits(p_mixE, "ggplot"))
ok("build_panelE() keeps BOTH the low-n and normal-n feature (nothing is dropped for low n)",
   setequal(p_mixE$data$feature, c("f_lo", "f_hi")))
ok("build_panelE() flags the low-carrier feature as low_n and leaves the normal one alone",
   isTRUE(p_mixE$data$low_n[p_mixE$data$feature == "f_lo"]) &&
   isFALSE(p_mixE$data$low_n[p_mixE$data$feature == "f_hi"]))

built_mixE <- try(ggplot2::ggplot_build(p_mixE), silent = TRUE)
ok("build_panelE() renders without error on a mix of low-n and normal-n rows",
   !inherits(built_mixE, "try-error"))
alpha_mixE <- built_mixE$data[[which(vapply(p_mixE$layers, function(l)
  inherits(l$geom, "GeomPoint"), logical(1)))]]$alpha
ok("the low-n point is rendered at reduced alpha and the normal-n point is unchanged (alpha 1)",
   isTRUE(all.equal(sort(alpha_mixE), sort(c(0.35, 1)))))

mix_forest <- tibble::tibble(
  analysis      = rep("A", 3),
  pair          = c("eQTL-pQTL", "eQTL-pQTL", "eQTL-sQTL"),
  feature       = c("f1", "f2", "f3"),
  odds_ratio    = c(2.0, 3.0, 1.5),
  conf.low      = c(1.1, 1.5, 0.9),
  conf.high     = c(3.6, 5.0, 2.4),
  p.value       = c(0.01, 0.001, 0.2),
  padj          = c(0.01, 0.001, 0.2),
  n_carriers    = c(6L, 25L, 40L),     # f1 is a lead, f2/f3 are normal
  flagged       = c(FALSE, FALSE, FALSE),
  estimate_type = "odds_ratio")

p_mix_forest <- enrichment_forest_panel(mix_forest, "A")
ok("enrichment_forest_panel() returns a ggplot object on a mix of low-n and normal-n rows",
   inherits(p_mix_forest, "ggplot"))

built_mix_forest <- try(ggplot2::ggplot_build(p_mix_forest), silent = TRUE)
ok("enrichment_forest_panel() renders without error on a mix of low-n and normal-n rows",
   !inherits(built_mix_forest, "try-error"))
alpha_mix_forest <- built_mix_forest$data[[which(vapply(p_mix_forest$layers, function(l)
  inherits(l$geom, "GeomPointrange"), logical(1)))]]$alpha
ok("the low-n row (f1, n_carriers = 6) is rendered at reduced alpha; the normal-n rows are alpha 1",
   isTRUE(all.equal(sort(alpha_mix_forest), sort(c(0.4, 1, 1)))))

cat("test_fig3_enrichment.R: all assertions passed\n")
