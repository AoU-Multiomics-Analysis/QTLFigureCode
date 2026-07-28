# Analysis B: continuous log MAF Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Switch Fig3's Analysis B from a binary `rare` (MAF < 1%) outcome to continuous log10 MAF, keeping the binary result as a supplementary column.

**Architecture:** `.enrich_pairs()` gains a per-analysis model family. Analysis A is unchanged (binomial, adjusted OR). Analysis B fits `lm(log_maf ~ feature + PIP_qtl + log10(cs_size_med))` and reports a fold-change in MAF. Both binary and continuous B results are returned so the supplement can show each.

**Tech Stack:** R, dplyr, purrr, stats. Tests via `tests/helpers.R`, discovered by `tests/run_all.R`.

## Global Constraints

- Run every script from the project root; scripts self-reference as lowercase `fig<n>/…`.
- Never `git add -A` — stage named paths only. `tools/build_notebooks.py` and the `*.ipynb` files carry the repo owner's unrelated work; do not stage them.
- Tests must pass with no controlled-access data: `Rscript tests/run_all.R`.
- **Analysis A is unchanged.** Panel E must render identically — it plots analysis A only.
- Analysis B's statistic is NOT an odds ratio and must never share an axis with A's.

---

### Task 1: Continuous log MAF for analysis B

**Files:**
- Modify: `Fig3/00_setup.R` (`.enrich_pairs()`, `fit_feature_enrichment()`)
- Modify: `Fig3/supp_coloc_enrichment.R` (axis label for the B facet)
- Modify: `Fig3/NOTES.md`
- Test: `tests/test_fig3_enrichment.R`

**Interfaces:**
- Consumes: `fit_feature_enrichment(d, outcome, features, min_carriers)`, `.enrich_pairs(v, min_carriers)`.
- Produces: `.enrich_pairs()` returns the existing columns plus `estimate_type` (`"odds_ratio"` or `"maf_fold_change"`). Analysis B rows carry `odds_ratio = 10^(-coef)` — the FOLD-CHANGE by which carriers' MAF is LOWER — with CIs transformed the same way, so a value > 1 means "carriers are rarer" and reads in the same direction as analysis A's enrichment.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_fig3_enrichment.R` before the final `cat()`:

```r
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript tests/test_fig3_enrichment.R`
Expected: FAIL — `estimate_type` does not exist.

- [ ] **Step 3: Implement**

In `fit_feature_enrichment()`, add a `family` argument (`"binomial"` or `"gaussian"`, default `"binomial"`).

For `"gaussian"`: fit `stats::lm(outcome ~ feature + PIP_qtl + log10(cs_size_med))`, take `confint()` on the coefficient, and transform estimate and both CI bounds by `10^(-x)` so that "carriers are rarer" yields a value ABOVE 1, matching analysis A's direction convention. Note the sign flip also SWAPS the CI bounds — after transforming, `conf.low` must still be the smaller number. Keep the carrier floor, the per-feature `tryCatch` isolation, and the named-drop messages exactly as they are.

In `.enrich_pairs()`, analysis B computes `log_maf = log10(max_af)` on rows with `max_af > 0`, calls `fit_feature_enrichment(..., family = "gaussian")`, and tags `estimate_type = "maf_fold_change"`. Analysis A tags `"odds_ratio"` and is otherwise untouched.

Guard the zero/negative case: `max_af == 0` gives `-Inf` and must be dropped before fitting, with the count reported in a message.

- [ ] **Step 4: Run tests**

Run: `Rscript tests/test_fig3_enrichment.R && Rscript tests/run_all.R`
Expected: both PASS.

- [ ] **Step 5: Update the supplement and NOTES.md**

In `Fig3/supp_coloc_enrichment.R`, the analysis-B facet's x-axis label becomes `"MAF fold change (carriers rarer →)"`; analysis A keeps `"Adjusted odds ratio"`. Both remain on a log10 axis with the reference line at 1.

Append to `Fig3/NOTES.md` under the panel E section: analysis B now uses continuous log MAF; record the measured justification (5 → 7 significant, 87% direction agreement, promoter/enhancer contrast crossing to significance, permutation-calibrated despite non-normal residuals) and the one loss (`non_coding_transcript_variant`, tail-specific).

- [ ] **Step 6: Verify against real data**

```bash
Rscript -e 'Sys.setenv(FIG3_INPUT="/Users/evinmpadhi/PaperFigures", FIG3_DATA=tempfile(), FIG3_OBJ=tempfile(), FIG3_OUTDIR=tempfile()); source("fig3/00_setup.R"); e <- coloc_feature_enrichment(refresh=TRUE); print(table(e$analysis, e$estimate_type)); print(subset(e, analysis=="B" & padj<0.05)[, c("pair","feature","odds_ratio","padj")])'
```

Expect analysis B to report `maf_fold_change` and to recover PLS, dELS, pELS among its significant features. Report the actual output.

- [ ] **Step 7: Commit**

```bash
git add Fig3/00_setup.R Fig3/supp_coloc_enrichment.R Fig3/NOTES.md tests/test_fig3_enrichment.R
git commit -m "Fig3: analysis B uses continuous log MAF instead of a binary rare label"
```
