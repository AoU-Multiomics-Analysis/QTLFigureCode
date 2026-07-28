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

cat("test_fig3_enrichment.R: all assertions passed\n")
