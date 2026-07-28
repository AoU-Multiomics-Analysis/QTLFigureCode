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

cat("test_fig3_enrichment.R: all assertions passed\n")
