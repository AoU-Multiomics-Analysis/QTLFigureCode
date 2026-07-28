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

# --- per-feature model fit --------------------------------------------------
set.seed(2026)
n <- 800
sim <- data.frame(
  PIP_qtl     = runif(n, 0.1, 0.9),
  cs_size_med = sample(2:40, n, replace = TRUE),
  # `hit` is strongly associated with the outcome; `dud` is not.
  hit         = rep(c(TRUE, FALSE), each = n / 2),
  dud         = sample(c(TRUE, FALSE), n, replace = TRUE),
  # `thin` has only 4 carriers -- below the floor of 10.
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
# `separated` clears the >=10-carrier floor (its smallest table(x, outcome)
# cell is well above 10) but is constructed so that INCLUDING it as a glm
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
# per-MAF-bin regressions, not for a 10-carrier floor applied within a single
# assay pair: min(table(feature, outcome)) can never reach 10 out of 9 total
# rows, so it is expected (not a bug) that no (analysis, pair) combination
# here survives fit_feature_enrichment()'s carrier floor. The assertions
# below therefore check the SHAPE of the result unconditionally, and check
# per-row properties only if any rows happen to exist -- exactly per the
# "do not depend on every combination succeeding" contract.
ok("the documented columns are all present",
   all(c("analysis", "pair", "feature", "odds_ratio", "conf.low", "conf.high",
         "p.value", "padj", "n_carriers", "flagged") %in% names(enr)))
ok("coloc_feature_enrichment() returns a well-formed tibble even when the fixture is too small for any feature to survive",
   is.data.frame(enr))
ok("only the documented analysis labels ever appear",
   all(unique(enr$analysis) %in% c("A", "B")))
ok("both analyses return the SAME columns",
   identical(sort(names(enr[enr$analysis == "A", ])),
             sort(names(enr[enr$analysis == "B", ]))))

if (nrow(enr)) {
  ok("both analyses are present", setequal(unique(enr$analysis), c("A", "B")))
  ok("padj is FDR within (analysis, pair), so never below the raw p",
     all(enr$padj >= enr$p.value - 1e-9))
  ok("flagged is logical", is.logical(enr$flagged))
  ok("odds ratios are finite -- the carrier floor rules out separation",
     all(is.finite(enr$odds_ratio)))
} else {
  message("  note: the fixture is too small for any (analysis, pair) ",
          "combination to clear the 10-carrier floor; coloc_feature_enrichment() ",
          "still returned a well-formed zero-row tibble of the documented shape.")
  ok("flagged column is logical-typed even with zero rows", is.logical(enr$flagged))
}

# --- .enrich_pairs() orchestration on a synthetic feature matrix ------------
# The real Fig3 fixture (9 rows) can never clear the 10-carrier floor, so the
# block above is schema-only -- deleting the FDR line, reverting `flagged`'s
# denominator, or dropping the `rare` exclusion from B would all leave `enr`
# an identical 0-row tibble there. This block builds a small, fully
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
  tibble::tibble(
    type = pair, Coloc = Coloc, hit = hit, dud1 = dud1, dud2 = dud2, dud3 = dud3,
    rare = rare,
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

cat("test_fig3_enrichment.R: all assertions passed\n")
