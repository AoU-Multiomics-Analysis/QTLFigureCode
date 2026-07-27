# Fig3 Panel E — Colocalization Feature Enrichment Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fill Fig3's reserved panel E slot with a feature-enrichment analysis of colocalizing versus QTL-layer-specific variants, plus a standalone supplement carrying the full result set.

**Architecture:** Two binomial-GLM analyses over one binary feature set, both adjusted for fine-mapping power (`PIP_qtl` + `log10(cs_size_med)`). A widened variant→feature map feeds a memoized builder that returns one tidy frame; two plotting scripts consume that same frame. All new logic lives behind pure seams so it is testable without controlled-access data.

**Tech Stack:** R, data.table, dplyr, broom, ggplot2, patchwork. Tests are plain scripts using `tests/helpers.R` (`ok()` / `expect_error()`), discovered by `tests/run_all.R` via the `^test_.*\.R$` pattern.

## Global Constraints

- Run every script from the **project root**. Scripts self-reference as lowercase `fig3/…`.
- **Never `git add -A`** — stage named paths only. Inputs and derived `.rds` are participant-derived.
- Tests must pass with **no controlled-access data**: `Rscript tests/run_all.R`.
- Fixtures resolve every path through the manifest (`manifest_path()` / `fixture_path()` in `tests/fixtures.R`) — never hardcode a layout.
- Any change that alters a **result** goes in `Fig3/NOTES.md`.
- Feature names match **VEP's own spelling** exactly (`3_prime_UTR_variant`, not `_utr_`).
- Count floor is **≥10 carriers in the smaller group**; dropped features are **named in a message**, never silently discarded.
- Cut-points: `PHYLOP_CONSERVED <- 2`, `TSS_PROXIMAL_BP <- 2000`.

---

## File Structure

| File | Responsibility |
|---|---|
| `common/vep.R` | **Create.** `expand_vep_consequences()` — one-hot expansion of comma-separated VEP terms. Promoted from Fig1 so both figures share one copy. |
| `Fig1/00_setup.R` | **Modify.** Delete the local `.expand_vep_consequences()`, source `common/vep.R` instead. Behaviour unchanged. |
| `Fig3/00_setup.R` | **Modify.** Widen the SuSiE select; `.variant_af_map()` → `.variant_feature_map()`; add feature constants, `build_feature_matrix()`, `fit_feature_enrichment()`, `coloc_feature_enrichment()`. |
| `Fig3/panels/panelE.R` | **Modify.** Replace the blank placeholder with the pair-contrast scatter. |
| `Fig3/supp_coloc_enrichment.R` | **Create.** Standalone full forest, both analyses, all three pairs. |
| `tests/fixtures.R` | **Modify.** Extend the Fig3 SuSiE fixture with annotation columns. |
| `tests/test_vep.R` | **Create.** Pins `expand_vep_consequences()` behaviour. |
| `tests/test_fig3_enrichment.R` | **Create.** Feature matrix, model fit, builder, flagging. |
| `Fig3/NOTES.md` | **Modify.** Record the filled slot and the open limitations. |

---

### Task 1: Promote the VEP expansion to `common/`

**Files:**
- Create: `common/vep.R`
- Modify: `Fig1/00_setup.R` (remove local copy, add `source()`)
- Test: `tests/test_vep.R`

**Interfaces:**
- Consumes: nothing.
- Produces: `expand_vep_consequences(d)` — takes a data frame with a `consequence` character column of comma-separated VEP terms; returns the same frame with one **logical** column added per distinct term. Note the **leading dot is dropped** on promotion (it is now a shared public helper, not a file-local one).

- [ ] **Step 1: Write the failing test**

Create `tests/test_vep.R`:

```r
# Behaviour of the shared VEP consequence expansion (common/vep.R).
# Run from the project root:  Rscript tests/test_vep.R
suppressPackageStartupMessages({ library(dplyr) })
source("common/vep.R")
source("tests/helpers.R")

d <- data.frame(
  variant     = paste0("v", 1:4),
  consequence = c("intron_variant,NMD_transcript_variant",
                  "missense_variant",
                  ".",
                  "splice_region_variant,intron_variant"),
  PLS         = c(TRUE, FALSE, TRUE, FALSE),
  stringsAsFactors = FALSE)
e <- expand_vep_consequences(d)

ok("adds one column per distinct term",
   all(c("intron_variant", "NMD_transcript_variant", "missense_variant",
         "splice_region_variant") %in% names(e)))
ok("added columns are logical",
   all(vapply(e[c("intron_variant", "missense_variant")], is.logical, logical(1))))
ok("a multi-term row sets EVERY one of its terms",
   isTRUE(e$splice_region_variant[4]) && isTRUE(e$intron_variant[4]))
ok("a term absent from a row is FALSE, not NA",
   isFALSE(e$missense_variant[1]))
ok("VEP's '.' placeholder produces no column", !"." %in% names(e))
ok("an existing column is never clobbered", identical(e$PLS, d$PLS))

# A frame with no consequence column passes straight through.
ok("no consequence column is a no-op",
   identical(expand_vep_consequences(data.frame(x = 1)), data.frame(x = 1)))
ok("zero rows is a no-op",
   nrow(expand_vep_consequences(d[0, ])) == 0L)

cat("test_vep.R: all assertions passed\n")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript tests/test_vep.R`
Expected: FAIL — `cannot open file 'common/vep.R'`

- [ ] **Step 3: Create `common/vep.R` and rewire Fig1**

Create `common/vep.R` by moving the function body currently in `Fig1/00_setup.R`, dropping the leading dot from the name:

```r
# =============================================================================
# VEP consequence expansion, shared by Fig1 panel H and Fig3 panel E.
#
# VEP reports every consequence for a variant as ONE comma-separated string
# ("intron_variant,NMD_transcript_variant"), and ~20% of rows carry more than
# one term. Enrichment code picks annotation columns up via where(is.logical),
# so unless each term becomes its own logical column the whole coding / splice
# / UTR / noncoding-RNA half of the annotation space is silently absent.
#
# One-hot expansion is done by matrix construction rather than separate_rows()
# + pivot_wider(): these tables carry ~80 columns, and pivot_wider would treat
# every one of them as an id column.
# =============================================================================
expand_vep_consequences <- function(d) {
  if (!"consequence" %in% names(d) || !nrow(d)) return(d)
  terms <- lapply(strsplit(as.character(d$consequence), ",", fixed = TRUE),
                  function(x) unique(trimws(x)))
  lvl <- setdiff(unique(unlist(terms)), c("", ".", NA))
  if (!length(lvl)) return(d)
  # Never clobber a column that already exists (the cCRE flags are already
  # one-hot); a duplicate name would make select() ambiguous downstream.
  lvl <- setdiff(lvl, names(d))
  if (!length(lvl)) return(d)
  m <- vapply(lvl, function(t) vapply(terms, function(x) t %in% x, logical(1)),
              logical(length(terms)))
  dplyr::bind_cols(d, as.data.frame(m, stringsAsFactors = FALSE))
}
```

In `Fig1/00_setup.R`: delete the whole `.expand_vep_consequences <- function(d) { … }` block **and its comment header**, add `source("common/vep.R")` immediately after the existing `source("common/manifest.R")` line, and change the single call site inside `ComputeFinemappingAnnotationEnrichment()` from `.expand_vep_consequences()` to `expand_vep_consequences()`.

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript tests/test_vep.R && Rscript tests/run_all.R`
Expected: both PASS. Fig1's panel H behaviour is unchanged — the function is identical, only its location and name changed.

- [ ] **Step 5: Commit**

```bash
git add common/vep.R Fig1/00_setup.R tests/test_vep.R
git commit -m "Promote the VEP consequence expansion to common/"
```

---

### Task 2: Widen the variant→feature map

**Files:**
- Modify: `Fig3/00_setup.R` (SuSiE select, `.variant_af_map()` → `.variant_feature_map()`, new constants)
- Modify: `tests/fixtures.R` (annotation columns on the Fig3 SuSiE fixture)
- Test: `tests/test_fig3_enrichment.R`

**Interfaces:**
- Consumes: `susie_finemapping()` (existing), `expand_vep_consequences()` from Task 1.
- Produces:
  - `FEATURE_CCRE_COLS` — character(7), the cCRE flag names.
  - `.variant_feature_map(fm)` — one row per variant, columns: `variant`, `AF_POP_COLS`, `FEATURE_CCRE_COLS`, `consequence`, `phylop`, `distTSS`, `gene_type`, `cs_size_med`. Retains the existing conflicting-AF error.

**Note on `cs_size_med`:** `cs_size` is a property of a credible set, and one variant can sit in several across genes and assays. The map is variant-level, so it carries the **median** `cs_size` over that variant's rows. This is a plan-level decision the spec left open; it keeps the map one-row-per-variant, which the many-to-one join guard in `coloc_variant_level()` depends on.

- [ ] **Step 1: Write the failing test**

Create `tests/test_fig3_enrichment.R`:

```r
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript tests/test_fig3_enrichment.R`
Expected: FAIL — `could not find function ".variant_feature_map"`

- [ ] **Step 3: Widen the fixture, the SuSiE select, and the map**

In `tests/fixtures.R`, inside `make_fig3_fixture()`, add annotation columns to all three SuSiE frames before they are written. Add this helper next to the existing `af_for()` and call it on `eqtl`, `sqtl` and `pqtl` (`eqtl <- cbind(eqtl, annot_for(eqtl$variant))`, etc.):

```r
  # Annotation columns for the enrichment feature set. Deterministic from the
  # variant index so assertions can target specific rows. Variant 1 is PLS +
  # missense + conserved + TSS-proximal; variant 2 carries two VEP terms.
  annot_for <- function(v) {
    i <- seq_along(v)
    data.frame(
      PLS             = i %% 2 == 1,
      pELS            = i %% 3 == 0,
      dELS            = i %% 4 == 0,
      `CTCF.bound`    = i %% 5 == 0,
      `CTCF.only`     = FALSE,
      `DNase.H3K4me3` = i %% 6 == 0,
      FANTOM5         = i %% 7 == 0,
      consequence     = ifelse(i == 1, "missense_variant",
                        ifelse(i == 2, "intron_variant,NMD_transcript_variant",
                               "intron_variant")),
      phylop          = ifelse(i %% 2 == 1, 3.5, 0.2),
      distTSS         = ifelse(i %% 2 == 1, 500L, 50000L),
      gene_type       = ifelse(i %% 3 == 0, "lncRNA", "protein_coding"),
      cs_size         = 2L + (i %% 4),
      check.names     = FALSE, stringsAsFactors = FALSE)
  }
```

In `Fig3/00_setup.R`, add the constants next to `AF_POP_COLS`:

```r
# cCRE flags carried through for the panel E enrichment. Already logical in the
# SuSiE tables.
FEATURE_CCRE_COLS <- c("PLS", "pELS", "dELS", "CTCF.bound", "CTCF.only",
                       "DNase.H3K4me3", "FANTOM5")
# Extra columns the enrichment needs beyond the base fine-mapping fields.
FEATURE_EXTRA_COLS <- c("consequence", "phylop", "distTSS", "gene_type", "cs_size")
```

Widen the three `select()` calls in `susie_finemapping()` to keep them — append `any_of(c(FEATURE_CCRE_COLS, FEATURE_EXTRA_COLS))` to each. `any_of()` (not `all_of()`) so a table lacking one optional annotation still loads.

Replace `.variant_af_map()` with:

```r
# One row per variant, carrying allele frequencies AND the annotation features
# panel E's enrichment needs.
#
# cs_size is a property of a CREDIBLE SET and one variant can sit in several
# across genes and assays, so the map carries the median over that variant's
# rows. Keeping exactly one row per variant is what the many-to-one join guard
# in coloc_variant_level() relies on.
.variant_feature_map <- function(fm) {
  keep <- c("variant", AF_POP_COLS, FEATURE_CCRE_COLS, FEATURE_EXTRA_COLS)
  raw <- bind_rows(lapply(list(fm$eQTL, fm$sQTL, fm$pQTL),
                          function(d) dplyr::select(d, any_of(keep)))) %>%
    # The eQTL reader already repairs "chrchr"; normalise all three here so a
    # stray prefix cannot split one variant into two rows.
    mutate(variant = str_replace(variant, "chrchr", "chr"),
           across(all_of(AF_POP_COLS), as.numeric))

  size <- raw %>%
    group_by(variant) %>%
    summarise(cs_size_med = stats::median(as.numeric(cs_size), na.rm = TRUE),
              .groups = "drop")

  af <- raw %>% select(-any_of("cs_size"))
  af <- af[rowSums(!is.na(af[AF_POP_COLS])) > 0, , drop = FALSE] %>% distinct()

  # A duplicate here means two SuSiE tables disagree about the same variant's
  # variant-level annotations, which is a genuine upstream inconsistency rather
  # than something to silently pick a winner for.
  dup <- unique(af$variant[duplicated(af$variant)])
  if (length(dup))
    stop(length(dup), " variant(s) carry CONFLICTING variant-level values ",
         "across the SuSiE fine-mapping tables -- e.g. '", dup[1],
         "'. These are properties of the variant alone; resolve the ",
         "disagreement upstream.", call. = FALSE)
  left_join(af, size, by = "variant")
}
```

Update the single call site in `coloc_variant_level()` from `.variant_af_map(susie_finemapping())` to `.variant_feature_map(susie_finemapping())`.

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript tests/test_fig3_enrichment.R && Rscript tests/run_all.R`
Expected: both PASS. The pre-existing "conflicting per-assay allele frequencies error the variant-level join" assertion in `tests/test_fig3_builders.R` still passes — the guard is retained, only its message widened.

- [ ] **Step 5: Commit**

```bash
git add Fig3/00_setup.R tests/fixtures.R tests/test_fig3_enrichment.R
git commit -m "Fig3: widen the variant map to carry enrichment features"
```

---

### Task 3: Build the binary feature matrix

**Files:**
- Modify: `Fig3/00_setup.R`
- Test: `tests/test_fig3_enrichment.R` (append)

**Interfaces:**
- Consumes: `expand_vep_consequences()` (Task 1), `.variant_feature_map()` (Task 2).
- Produces: `build_feature_matrix(d)` — takes a frame with `consequence`, `phylop`, `distTSS`, `gene_type`, `max_af` and the cCRE flags; returns it with four binarized logical columns added (`rare`, `conserved`, `tss_proximal`, `protein_coding`) and VEP terms expanded. Also `feature_names(d)` — the character vector of logical column names eligible as features.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_fig3_enrichment.R`, before the final `cat()`:

```r
# --- feature matrix ---------------------------------------------------------
fmx <- build_feature_matrix(data.frame(
  variant     = paste0("v", 1:3),
  consequence = c("missense_variant", "intron_variant", "."),
  phylop      = c(3.5, 0.2, 1.0),
  distTSS     = c(500, 50000, -100),
  gene_type   = c("protein_coding", "lncRNA", "protein_coding"),
  max_af      = c(0.001, 0.20, 0.05),
  PLS         = c(TRUE, FALSE, TRUE),
  stringsAsFactors = FALSE))

ok("rare uses MAF < 0.01",        identical(fmx$rare, c(TRUE, FALSE, FALSE)))
ok("conserved uses phyloP > 2",   identical(fmx$conserved, c(TRUE, FALSE, FALSE)))
ok("tss_proximal uses |distTSS| < 2000 and is sign-agnostic",
   identical(fmx$tss_proximal, c(TRUE, FALSE, TRUE)))
ok("protein_coding reads gene_type",
   identical(fmx$protein_coding, c(TRUE, FALSE, TRUE)))
ok("VEP terms are expanded into the matrix", "missense_variant" %in% names(fmx))
ok("feature_names() returns only logical columns",
   all(vapply(fmx[feature_names(fmx)], is.logical, logical(1))))
ok("feature_names() excludes the outcome-ish columns",
   !any(c("variant", "consequence", "gene_type") %in% feature_names(fmx)))
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript tests/test_fig3_enrichment.R`
Expected: FAIL — `could not find function "build_feature_matrix"`

- [ ] **Step 3: Implement**

Add to `Fig3/00_setup.R`:

```r
# Cut-points for the binarized covariates. Conventional defaults, NOT tuned --
# if a result turns on the exact value that is a finding about fragility and
# belongs in NOTES.md.
PHYLOP_CONSERVED <- 2
TSS_PROXIMAL_BP  <- 2000

# Columns that are logical but are NOT features (outcomes or bookkeeping).
NON_FEATURE_LOGICAL <- c("Coloc", "coloc", "flagged")

# Turn the raw annotation columns into the binary feature set: four binarized
# covariates plus one logical column per VEP consequence term.
build_feature_matrix <- function(d) {
  d %>%
    mutate(
      rare           = !is.na(max_af)    & max_af   <  0.01,
      conserved      = !is.na(phylop)    & phylop   >  PHYLOP_CONSERVED,
      # Sign-agnostic: distTSS is signed by strand, proximity is not.
      tss_proximal   = !is.na(distTSS)   & abs(distTSS) < TSS_PROXIMAL_BP,
      protein_coding = !is.na(gene_type) & gene_type == "protein_coding") %>%
    expand_vep_consequences()
}

# Logical columns eligible as enrichment features.
feature_names <- function(d)
  setdiff(names(d)[vapply(d, is.logical, logical(1))], NON_FEATURE_LOGICAL)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript tests/test_fig3_enrichment.R && Rscript tests/run_all.R`
Expected: both PASS.

- [ ] **Step 5: Commit**

```bash
git add Fig3/00_setup.R tests/test_fig3_enrichment.R
git commit -m "Fig3: binary feature matrix for the coloc enrichment"
```

---

### Task 4: The per-feature model fit (pure seam)

**Files:**
- Modify: `Fig3/00_setup.R`
- Test: `tests/test_fig3_enrichment.R` (append)

**Interfaces:**
- Consumes: nothing from earlier tasks — takes a plain frame, so it is testable in isolation.
- Produces: `fit_feature_enrichment(d, outcome, features, min_carriers = 10)` — `d` has a logical `outcome` column, numeric `PIP_qtl` and `cs_size_med`, and the logical `features`. Returns a tibble with one row per **surviving** feature: `feature`, `odds_ratio`, `conf.low`, `conf.high`, `p.value`, `n_carriers`. Features below the floor are dropped and **named in a message**.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_fig3_enrichment.R`:

```r
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript tests/test_fig3_enrichment.R`
Expected: FAIL — `could not find function "fit_feature_enrichment"`

- [ ] **Step 3: Implement**

Add to `Fig3/00_setup.R`:

```r
# Adjusted enrichment of each feature on a binary outcome.
#
# The PIP + log10(cs_size) adjustment is the point of the design: colocalization
# requires good fine-mapping in BOTH layers, so an unadjusted test conflates
# biology with detectability.
#
# A feature is fitted only with at least `min_carriers` observations in the
# SMALLER of its two groups. Below that, glm returns separated or wildly
# unstable estimates -- the infinite odds ratios that silently vanished from
# Fig1 panel H. Dropped features are named in a message rather than disappearing.
fit_feature_enrichment <- function(d, outcome, features, min_carriers = 10) {
  empty <- tibble::tibble(feature = character(), odds_ratio = numeric(),
                          conf.low = numeric(), conf.high = numeric(),
                          p.value = numeric(), n_carriers = integer())
  y <- as.logical(d[[outcome]])
  keep <- character(0); dropped <- character(0)
  for (f in features) {
    x <- as.logical(d[[f]])
    tab <- table(x, y)
    if (nrow(tab) == 2L && ncol(tab) == 2L && min(tab) >= min_carriers)
      keep <- c(keep, f) else dropped <- c(dropped, f)
  }
  if (length(dropped))
    message("  ", length(dropped), " feature(s) below the ", min_carriers,
            "-carrier floor, not tested: ", paste(dropped, collapse = ", "))
  if (!length(keep)) return(empty)

  purrr::map_dfr(keep, function(f) {
    form <- stats::as.formula(
      sprintf("%s ~ `%s` + PIP_qtl + log10(cs_size_med)", outcome, f))
    fit <- stats::glm(form, data = d, family = stats::binomial())
    ci  <- suppressMessages(stats::confint(fit))
    co  <- summary(fit)$coefficients
    row <- paste0("`", f, "`TRUE")
    if (!row %in% rownames(co)) row <- grep(f, rownames(co), fixed = TRUE, value = TRUE)[1]
    tibble::tibble(
      feature    = f,
      odds_ratio = unname(exp(stats::coef(fit)[row])),
      conf.low   = unname(exp(ci[row, 1])),
      conf.high  = unname(exp(ci[row, 2])),
      p.value    = unname(co[row, 4]),
      n_carriers = sum(as.logical(d[[f]]), na.rm = TRUE))
  })
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript tests/test_fig3_enrichment.R && Rscript tests/run_all.R`
Expected: both PASS.

- [ ] **Step 5: Commit**

```bash
git add Fig3/00_setup.R tests/test_fig3_enrichment.R
git commit -m "Fig3: adjusted per-feature enrichment fit with a carrier floor"
```

---

### Task 5: The memoized builder across analyses and pairs

**Files:**
- Modify: `Fig3/00_setup.R`
- Test: `tests/test_fig3_enrichment.R` (append)

**Interfaces:**
- Consumes: `coloc_variant_level()`, `build_feature_matrix()`, `feature_names()`, `fit_feature_enrichment()`.
- Produces: `coloc_feature_enrichment(refresh = FALSE)` — memoized via `build_or_cache("coloc_feature_enrichment", …)`. Returns a tibble: `analysis` (`"A"`/`"B"`), `pair`, `feature`, `odds_ratio`, `conf.low`, `conf.high`, `p.value`, `padj`, `n_carriers`, `flagged`.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_fig3_enrichment.R`:

```r
# --- builder ----------------------------------------------------------------
enr <- coloc_feature_enrichment(refresh = TRUE)

ok("both analyses are present", setequal(unique(enr$analysis), c("A", "B")))
ok("the documented columns are all present",
   all(c("analysis", "pair", "feature", "odds_ratio", "conf.low", "conf.high",
         "p.value", "padj", "n_carriers", "flagged") %in% names(enr)))
ok("both analyses return the SAME columns",
   identical(sort(names(enr[enr$analysis == "A", ])),
             sort(names(enr[enr$analysis == "B", ]))))
ok("padj is FDR within (analysis, pair), so never below the raw p",
   all(enr$padj >= enr$p.value - 1e-9))
ok("flagged is logical", is.logical(enr$flagged))
ok("odds ratios are finite -- the carrier floor rules out separation",
   all(is.finite(enr$odds_ratio)))
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript tests/test_fig3_enrichment.R`
Expected: FAIL — `could not find function "coloc_feature_enrichment"`

- [ ] **Step 3: Implement**

Add to `Fig3/00_setup.R`:

```r
# A pair is flagged underpowered for an analysis when fewer than this fraction
# of its candidate features clear the carrier floor.
UNDERPOWERED_FRAC <- 0.5

# Adjusted feature enrichment for both analyses, per assay pair.
#
#   A  what makes a variant shared    coloc vs non-coloc, within pair
#   B  what makes a shared variant rare   rare vs common, within Coloc == TRUE
#
# Analysis B is attempted for every pair; where a pair has too few rare
# colocalizing variants the carrier floor empties it and the underpowered flag
# fires on its own. The criterion decides -- nothing is excluded a priori.
coloc_feature_enrichment <- function(refresh = FALSE)
  build_or_cache("coloc_feature_enrichment", refresh = refresh, {
    vl <- coloc_variant_level()
    # Fold ALT-allele frequencies to MAF, then take the per-variant maximum
    # across ancestries. Two steps, matching coloc_model_input() exactly --
    # inlining across() into do.call(pmax, ...) does not evaluate the way it
    # reads.
    folded <- vl %>%
      select(all_of(AF_POP_COLS)) %>%
      mutate(across(everything(), ~ ifelse(. > .5, 1 - ., .)))
    v <- vl %>%
      mutate(max_af = do.call(pmax, as.list(folded))) %>%
      build_feature_matrix()

    feats <- feature_names(v)
    message("  ", length(feats), " candidate features")

    specs <- tidyr::expand_grid(analysis = c("A", "B"), pair = unique(v$type))

    purrr::pmap_dfr(specs, function(analysis, pair) {
      d <- v %>% filter(type == pair)
      if (analysis == "A") {
        d$.outcome <- as.logical(d$Coloc)
      } else {
        d <- d %>% filter(Coloc == TRUE)
        d$.outcome <- as.logical(d$rare)
      }
      message("  analysis ", analysis, " / ", pair, " (n = ", nrow(d), ")")
      res <- fit_feature_enrichment(d, ".outcome", setdiff(feats, "rare"))
      # Coverage assertion: a shrinking feature set must not look like a clean run.
      n_expected <- length(setdiff(feats, "rare"))
      if (nrow(res) < n_expected)
        message("    ", nrow(res), "/", n_expected, " features tested")
      if (!nrow(res)) return(NULL)
      res %>%
        mutate(analysis = analysis, pair = pair,
               padj = stats::p.adjust(p.value, method = "fdr"),
               flagged = nrow(res) < UNDERPOWERED_FRAC * n_expected)
    })
  })
```

Note `rare` is excluded from the feature list — it is a *feature* in analysis A but the *outcome* in B, and regressing a variable on itself is meaningless. Excluding it from both keeps one feature list and avoids the trap.

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript tests/test_fig3_enrichment.R && Rscript tests/run_all.R`
Expected: both PASS.

- [ ] **Step 5: Commit**

```bash
git add Fig3/00_setup.R tests/test_fig3_enrichment.R
git commit -m "Fig3: memoized coloc feature-enrichment builder"
```

---

### Task 6: Panel E — the pair-contrast scatter

**Files:**
- Modify: `Fig3/panels/panelE.R`
- Modify: `Fig3/NOTES.md`
- Test: `tests/test_fig3_enrichment.R` (append)

**Interfaces:**
- Consumes: `coloc_feature_enrichment()`.
- Produces: `Fig3/objects/E.rds` — a ggplot. `assemble.R` is untouched.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_fig3_enrichment.R`:

```r
# --- panel E ----------------------------------------------------------------
pE <- build_panelE(enr)
ok("panel E is a ggplot", inherits(pE, "ggplot"))
ok("panel E plots analysis A only", all(pE$data$analysis == "A"))
ok("panel E is a two-pair contrast, one row per feature",
   !any(duplicated(pE$data$feature)))
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript tests/test_fig3_enrichment.R`
Expected: FAIL — `could not find function "build_panelE"`

- [ ] **Step 3: Implement**

Add `build_panelE()` to `Fig3/00_setup.R` (kept beside the other builders so the panel script stays thin):

```r
# Pair-contrast scatter: adjusted log-OR in eQTL-pQTL against eQTL-sQTL.
# Points on the diagonal are generic shared-signal markers; off-diagonal ones
# distinguish protein- from splicing-colocalization. sQTL-pQTL is structurally
# absent -- this plot has two axes -- and appears in the supplement instead.
build_panelE <- function(enr) {
  wide <- enr %>%
    filter(analysis == "A", pair %in% c("eQTL-pQTL", "eQTL-sQTL")) %>%
    mutate(log_or = log2(odds_ratio), sig = padj < 0.05) %>%
    select(analysis, feature, pair, log_or, sig) %>%
    tidyr::pivot_wider(names_from = pair,
                       values_from = c(log_or, sig)) %>%
    filter(!is.na(`log_or_eQTL-pQTL`), !is.na(`log_or_eQTL-sQTL`))

  ggplot(wide, aes(`log_or_eQTL-pQTL`, `log_or_eQTL-sQTL`)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey85") +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey85") +
    geom_point(aes(fill = `sig_eQTL-pQTL` & `sig_eQTL-sQTL`),
               shape = 21, size = 3.4, stroke = 0.4, colour = "grey25") +
    ggrepel::geom_text_repel(aes(label = feature), size = 3, max.overlaps = 12,
                             segment.colour = "grey70") +
    scale_fill_manual(values = c(`TRUE` = "#364B63", `FALSE` = "grey85"),
                      name = "FDR < 0.05 in both") +
    labs(x = "Adjusted log2 OR — eQTL-pQTL",
         y = "Adjusted log2 OR — eQTL-sQTL") +
    theme(legend.position = "bottom")
}
```

Replace the body of `Fig3/panels/panelE.R`:

```r
# Panel E — feature enrichment of colocalizing vs layer-specific variants,
# contrasted between assay pairs. See docs/superpowers/specs/
# 2026-07-27-fig3-panelE-coloc-enrichment-design.md
# Inputs: enloc_variant, susie_eqtl, susie_sqtl, susie_pqtl
source("fig3/00_setup.R")

pE <- build_panelE(coloc_feature_enrichment())

saveRDS(pE, file.path(OBJ_DIR, "E.rds"))
message("Panel E built.")
```

Append to `Fig3/NOTES.md`, replacing the `## Panel E — reserved slot (undecided)` section:

```markdown
## Panel E — colocalization feature enrichment

**This fills a previously blank slot, so Figure 3 gains a panel.**

Panel E contrasts the adjusted enrichment of binary variant features between
assay pairs: x = eQTL-pQTL, y = eQTL-sQTL, both as log2 odds ratios from
`coloc ~ feature + PIP_qtl + log10(cs_size_med)`. Design and rationale:
`docs/superpowers/specs/2026-07-27-fig3-panelE-coloc-enrichment-design.md`.

Limitations carried deliberately:

- **MAF is not adjusted for.** Rare-vs-common is a tested feature, so adjusting
  for it would condition away the effect being measured.
- **`cs_size_med` is a per-variant median.** `cs_size` belongs to a credible
  set and a variant can sit in several; the median keeps the feature map at one
  row per variant, which `coloc_variant_level()`'s join guard requires.
- **`max_af` NA behaviour is inherited** from `coloc_model_input()` — a variant
  missing any one per-ancestry frequency drops out.
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript tests/test_fig3_enrichment.R && Rscript tests/run_all.R`
Expected: both PASS, including the existing `test_fig3_assembly.R` tag test — `assemble.R` is unchanged and E is still one object.

- [ ] **Step 5: Commit**

```bash
git add Fig3/00_setup.R Fig3/panels/panelE.R Fig3/NOTES.md tests/test_fig3_enrichment.R
git commit -m "Fig3: fill panel E with the coloc enrichment contrast"
```

---

### Task 7: The supplement, and the notebook

**Files:**
- Create: `Fig3/supp_coloc_enrichment.R`
- Modify: `tools/build_notebooks.py`, `Fig3/Fig3.ipynb`
- Test: manual run (the supplement is standalone and not covered by `run_all.R`, matching `supp_coloc_qc.R`)

**Interfaces:**
- Consumes: `coloc_feature_enrichment()`.
- Produces: `FigS_ColocEnrichment.pdf` / `.png` in `OUTDIR`.

- [ ] **Step 1: Write the supplement**

Create `Fig3/supp_coloc_enrichment.R`:

```r
# =============================================================================
# Supplemental — full colocalization feature enrichment. NOT part of Figure 3.
# Produces FigS_ColocEnrichment.pdf:
#   analysis A — what makes a variant shared (coloc vs non-coloc, within pair)
#   analysis B — what makes a shared variant rare (within Coloc == TRUE)
#
# Panel E of the main figure shows only analysis A, and only the two
# well-powered pairs. Everything else lives here.
#
# Run from the project root:  Rscript fig3/supp_coloc_enrichment.R
# =============================================================================
source("fig3/00_setup.R")

enr <- coloc_feature_enrichment()

ANALYSIS_LABELS <- c(A = "A — shared vs layer-specific",
                     B = "B — rare vs common, among shared")

forest <- function(d, title) {
  ggplot(d, aes(x = odds_ratio, y = reorder(feature, odds_ratio),
                colour = flagged)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
    geom_pointrange(aes(xmin = conf.low, xmax = conf.high), size = 0.35) +
    facet_wrap(~pair, nrow = 1) +
    scale_x_continuous(trans = "log10") +
    # Grey = the pair was flagged underpowered for this analysis.
    scale_colour_manual(values = c(`FALSE` = "#364B63", `TRUE` = "grey70"),
                        guide = "none") +
    labs(x = "Adjusted odds ratio", y = NULL, title = title) +
    theme(plot.title = element_text(size = SZ_STRIP, face = "bold"))
}

figS <- forest(dplyr::filter(enr, analysis == "A"), ANALYSIS_LABELS[["A"]]) /
        forest(dplyr::filter(enr, analysis == "B"), ANALYSIS_LABELS[["B"]])

ggsave(file.path(OUTDIR, "FigS_ColocEnrichment.pdf"), figS,
       width = 14, height = 12, dpi = DPI, limitsize = FALSE)
ggsave(file.path(OUTDIR, "FigS_ColocEnrichment.png"), figS,
       width = 14, height = 12, dpi = DPI, limitsize = FALSE)
message("Assembled: ", file.path(OUTDIR, "FigS_ColocEnrichment.pdf"))
```

- [ ] **Step 2: Verify the supplement runs**

Run:
```bash
FIG3_INPUT=/Users/evinmpadhi/PaperFigures Rscript fig3/supp_coloc_enrichment.R
```
Expected: writes `FigS_ColocEnrichment.pdf` and `.png`. Confirm the file exists and is non-empty before claiming success.

- [ ] **Step 3: Register the supplement in the notebook generator**

In `tools/build_notebooks.py`, in the Fig3 config's `steps` list, add after the existing `supp` entry:

```python
     ("supp","Fig3/supp_coloc_enrichment.R",None,"figS",
        "**Supplemental (colocalization feature enrichment).** Adjusted odds ratios "
        "for every binary variant feature, both analyses, all three assay pairs. "
        "Panel E of the main figure carries only analysis A and the two "
        "well-powered pairs."),
```

- [ ] **Step 4: Regenerate and verify the notebook round-trips**

Run:
```bash
python3 tools/build_notebooks.py Fig3
a=$(shasum Fig3/Fig3.ipynb | cut -d' ' -f1)
python3 tools/build_notebooks.py Fig3
b=$(shasum Fig3/Fig3.ipynb | cut -d' ' -f1)
[ "$a" = "$b" ] && echo STABLE || echo NOT-STABLE
Rscript tests/run_all.R
```
Expected: `STABLE`, and the full suite passes.

- [ ] **Step 5: Commit**

```bash
git add Fig3/supp_coloc_enrichment.R tools/build_notebooks.py Fig3/Fig3.ipynb
git commit -m "Fig3: supplemental full coloc enrichment forest"
```

---

## Verification

Before claiming the plan complete, run and confirm output:

```bash
Rscript tests/run_all.R                    # every test file passes
FIG3_INPUT=/Users/evinmpadhi/PaperFigures Rscript fig3/run_all.R
FIG3_INPUT=/Users/evinmpadhi/PaperFigures Rscript fig3/supp_coloc_enrichment.R
git status --short                          # no participant-derived file staged
```

Panel E must render as a populated scatter, not the previous blank. `Fig3/Figure3.pdf` still carries exactly 5 tags.
