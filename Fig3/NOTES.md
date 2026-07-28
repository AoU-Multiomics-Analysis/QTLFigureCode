# Figure 3 — deferred work, open questions, and behaviour changes

## Behaviour changes from the original (each changes a scientific result)

These were found while rewriting the pipeline. None is a stylistic change —
each can move which rows are labelled colocalized, which variants enter a
panel, or what a rebuild actually reflects. Recorded here for whoever compares
against the pre-cleanup figure.

### 1. The colocalization-significance (FDR) threshold is now selected by rank, not reconstructed by subtraction

`add_fdr_threshold_label()` (`00_setup.R`) controls the Bayesian FDR at 5% by
sorting rows on their posterior probability, taking the local FDR `1 - p` at
each rank, and keeping the largest prefix whose running mean stays at or below
the target.

The original then rebuilt the cutoff as `1 - (1 - p)` — a round trip through
subtraction that is inexact in IEEE-754 floating point for roughly 18% of
realistic decimal probabilities below 0.5. That imprecision could make the row
that *defines* the FDR boundary fail its own `>=` cutoff test and be silently
excluded from the significant set it produced.

The rewrite instead reads the cutoff back as an actual observed element of the
probability column (`probs[ord[k]]`), never a reconstructed value. Two
consequences follow:

- Rows exactly tied at the cutoff are now **all** labelled the same way,
  rather than the outcome depending on sort order among ties.
- A row with a missing posterior is labelled `FALSE` rather than `NA`, so it
  can no longer poison a downstream `sum()`/`mean()` over the label column.

**Net effect: which rows are called colocalized can differ slightly from the
original at the FDR boundary**, in both the gene-, credible-set-, and
variant-level tables (`coloc_gene_level()`, `coloc_cs_level()`,
`coloc_variant_level()` all call this function).

### 2. The FDR cutoff itself carries a 1e-8 floating-point tolerance, which can admit one extra boundary row

`add_fdr_threshold_label()`'s `valid <- which(FDR <= control_level + 1e-8)` has
no counterpart in the original's `which(FDR <= control_level)`. The tolerance
guards against a specific IEEE-754 edge case: a running-mean FDR that is
mathematically exactly at `control_level` can round to a double a hair above
it, which would make `<=` reject a rank that should pass. Without the
tolerance, that rank (and, because `k <- max(valid)`, every rank has already
passed the same test up to it) is dropped from the significant prefix.

This is separate from item 1's read-back/tie-handling fix above: it changes
which *rank* qualifies as the cutoff in the first place, not how the cutoff
value is read back once chosen. **Net effect: at most one additional row,
whose running-mean FDR sits within 1e-8 above `control_level`, can be labelled
colocalized that the original's exact `<=` would have excluded.**

### 3. Panel D now excludes variants with no allele frequency, instead of silently binning them as "Rare"

`panels/panelD.R`'s original had `filter(!is.na(maf_bin))` — dead code,
because the `case_when()` that built `maf_bin` had a catch-all `TRUE ~
'Rare'` arm. A variant with an `NA` maximum allele frequency (`max_af`) never
matches `>= .05` or `>= .01`, falls through to the catch-all, and lands in the
**Rare** bin instead of being dropped. The `!is.na(maf_bin)` filter downstream
can never be `TRUE` for such a row and so removes nothing.

The rewrite filters on `!is.na(max_af)` **before** binning, so those variants
are excluded from the panel entirely rather than miscounted as rare. This
changes the Rare facet's N, its fitted slope, and — since the shared
`coord_cartesian()` limits are the 97th percentile of `|posterior_mean|`
pooled across all three facets — potentially the axis limits shown for all
three facets.

> **Alternative considered for `max_af` itself.** `coloc_model_input()`
> (`00_setup.R`) computes `max_af` with
> `do.call(pmax, as.list(folded))` — matching the original's
> `pmax(gvs_afr_af, ..., gvs_oth_af)` exactly: a variant missing even one of
> the seven per-ancestry frequencies gets `max_af = NA`, and is then excluded
> here by the filter above. Passing `na.rm = TRUE` to that `pmax` (i.e.
> changing the single `mutate(max_af = ...)` line in `coloc_model_input()`)
> is a defensible alternative: it would take the maximum over whichever
> ancestries *do* have a frequency, rather than discarding the variant
> outright, on the reasoning that a partial observation is still evidence
> about the variant's frequency. That is a deliberate scientific choice a
> researcher may prefer — it has not been made here, so `max_af` and every
> downstream quantity that reads it (`rare`, `log_maf`, panel D's bins, and
> all three of panel C's interaction estimates) currently reproduce the
> original's semantics exactly.

### 4. The credible-set colocalization graph now unions all three assay pairs, not just sQTL-pQTL — panel B is expected to change

The original's `Fig3.R` had a duplicated `bind_rows()` argument that listed
one assay pair twice and never bound eQTL-sQTL, so the graph behind panel B
and its inset (`coloc_cs_components()` / `coloc_cs_edges()`) was built from
sQTL-pQTL edges alone (see the README's [Provenance](README.md#provenance)
section and the regression test `tests/test_fig3_edges.R`). The rewrite fixes
that duplicate argument, so the graph is now the union of edges from all
three pairs.

**This is the single largest expected change to a published panel.** Every
count panel B derives from `coloc_cs_components()` — the distribution of
genes by number of colocalized signals, and the inset's per-assay
credible-set colocalization rate — will differ from the pre-cleanup figure.
That difference is the fix, not a regression.

### 5. New hard-stop cardinality guards can make the pipeline refuse to build on real data, with no override

The annotation join in `coloc_variant_level()` and both posterior joins in
`coloc_model_input()`'s `join_posteriors()` now pass
`relationship = "many-to-one"` to `left_join()`. The original performed the
same joins with no such assertion: a duplicated key on the right-hand side
(a repeated `SNP` in the annotation table, or a variant sitting in two
credible sets for the same gene) would silently fan the join out, inflating
row counts and, for the posterior joins, double-counting that variant in the
regression's sample size — the original never surfaced this even as a
warning.

The guards are the right fix, but the consequence is asymmetric: they convert
a silent correctness bug into a **hard `stop()`**, and there is no option to
proceed anyway. This pipeline has, to date, only ever run against the mock
fixtures in `tests/`, which are constructed to satisfy the guards by
assumption. On real *All of Us* data — which nothing upstream of these joins
validates — a duplicated key that the original would have silently mishandled
will now stop the pipeline from producing the figure at all. The researcher
should expect that possibility on the first real run, and should read the
`stop()` message (which names the offending table and, where findable, an
example offending key) rather than treat it as an unrelated crash.

### 6. `refresh = TRUE` does not cascade through the derived-data cache

Each builder in `00_setup.R` (`susie_finemapping()`, `coloc_gene_level()`,
`coloc_cs_level()`, `coloc_variant_level()`, `coloc_cs_components()`,
`coloc_model_input()`) memoizes under its own name in `fig3/data/` and calls
its dependencies as ordinary R calls. Passing `refresh = TRUE` to one builder
forces *that* builder to recompute, but any dependency it calls (e.g.
`coloc_model_input()` calling `coloc_variant_level()` and
`susie_finemapping()`) is called with its own default `refresh = FALSE` and so
still returns its cached, possibly stale, result.

There is no cascading invalidation. To rebuild everything after changing an
input file or a builder's logic, clear the whole cache directory
(`rm -rf fig3/data`) rather than passing `refresh = TRUE` to the one builder
you touched. (Also noted in the README's cache section.)

### 7. Allele frequencies now come from the SuSiE tables, not the AoU annotation table — panels C and D are expected to move

The `annotations` input
(`annotated/AoU.COMB.AC10.AN95.biallelic.v9.annotations.tsv.bgz`) has been
removed from `inputs.tsv` entirely. It was the largest input in the repo and
Fig3 used it for exactly one thing: the seven per-ancestry allele-frequency
columns (`AF_POP_COLS`), folded to MAF and reduced with `pmax()` to `max_af`
(and hence `rare` and `log_maf`). `chrom`/`pos`/`ref`/`alt` only built the
join key; nothing else in the table was ever read.

Those same seven columns ride along on the annotated SuSiE outputs, so
`coloc_variant_level()` (`00_setup.R`, previously the `fread(input_path(
"annotations"))` + `left_join(..., by = "SNP")` block) now sources them from
`susie_finemapping()` via the new `.variant_af_map()`, which pools the three
assay tables and reduces them to one row per variant.

**Why this can move a result — two distinct mechanisms:**

1. **Coverage.** The annotation table was genome-wide; the SuSiE tables are
   not. A variant in the FastENLOC `.snp.out` output that appears in no SuSiE
   table now gets no allele frequency, and is dropped downstream by the
   `max_af` NA filters (panel D's `filter(!is.na(max_af))`, panel C's NA
   filter). `coloc_variant_level()` emits a `note:` line reporting how many
   rows this affects — **check it on the first real rebuild.** Under the old
   code these variants were annotated and retained.

2. **`log_maf` is standardised over the pooled set.** In
   `coloc_model_input()`, `log_maf = as.numeric(scale(-log10(max_af)))` is
   computed *before* `filter(Coloc == TRUE)`, across all three pairs at once.
   `scale()`'s mean and SD therefore depend on exactly which variants carry a
   non-NA `max_af`. Any coverage change shifts the centring for **every**
   variant, so panel C's interaction coefficients and confidence intervals
   move even for variants whose own frequency is unchanged. This is why the
   change is recorded here rather than treated as a refactor.

Where both sources have a variant, the frequencies themselves are expected to
be identical — the SuSiE annotation step draws from the same AoU callset. The
folding is unchanged and still uses the seven per-ancestry columns, **not**
the precomputed `gvs_max_af`, which is an unfolded max ALT frequency and would
give a different answer (the fixture pins this: a variant with
`gvs_eur_af = 0.90` must fold to 0.10, leaving `max_af = 0.30`).

Because allele frequency is a property of the variant rather than of an assay,
`.variant_af_map()` hard-stops if two SuSiE tables disagree about the same
variant, rather than silently picking a winner — a disagreement would give one
variant two AF rows and fan out the variant-level join.

**Existing caches are stale.** Per note 6 there is no cascading invalidation,
so `rm -rf fig3/data` before the first rebuild or `coloc_variant_level()` will
return a cached frame still built from the annotation table.

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

### The carrier floor is 5, not 10

`fit_feature_enrichment()`'s `min_carriers` (and `.enrich_pairs()`'s own
default, which threads it through) was **10**, chosen against an earlier
figure that produced infinite odds ratios from separated fits. Measured
empirically on the real staged data at both values:

| | floor 10 | floor 5 |
|---|---|---|
| tests fitted | 106 | 124 |
| FDR-significant | 34 | **40** |
| features lost relative to floor 10 | — | **none** |

Every feature significant at floor 10 is still significant at floor 5, so
lowering the floor is a strict gain: six additional FDR-significant features,
with no loss. Four of the six are substantial (n = 18–44 carriers) and
include the largest effects in the study — `stop_gained` (OR 6.80),
`frameshift_variant` (OR 6.58), `splice_acceptor_variant` (OR 7.31). The
remaining two rest on very few carriers and are fragile:

- analysis B's `stop_gained`, n = 7, CI [1.71, 13.55]
- analysis B's `frameshift_variant`, n = 8, CI [2.68, 18.55]

Both are still worth showing — the effect sizes are large and consistent
with the truncating-variant biology of the four substantial features above —
but with only 7–8 carriers behind them they are **leads, not stable
estimates**. Rather than hide that distinction in the data (or drop the
floor further, or leave it at 10 and lose real findings), any row with
`n_carriers < 10` is rendered visually weaker (reduced alpha, `low_n` in
`build_panelE()`) in both panel E's scatter and the
`enrichment_forest_panel()` supplement, with a caption noting the
distinction, so a reader sees at a glance which points are estimates and
which are leads. The 10-carrier visual cutoff is deliberately kept at the
OLD statistical floor value, now repurposed as a legibility threshold rather
than an inclusion threshold.

### Analysis B now uses continuous log10 MAF instead of a binary `rare` outcome

Panel E's underlying enrichment fits two analyses per assay pair (see
`.enrich_pairs()`, `00_setup.R`): analysis A ("what makes a variant shared",
`Coloc ~ feature + PIP_qtl + log10(cs_size_med)`) feeds panel E and is
**unchanged by this note**. Analysis B ("what distinguishes shared variants
that are rare") originally dichotomised minor allele frequency at 1%
(`rare = max_af < .01`) and fit the same binomial model on that binary label.

Measuring the real staged data showed that cut is badly placed: the median
colocalizing variant sits at MAF ≈ 20%, and only 6–7% of colocalizing
variants fall below the 1% threshold analysis B was splitting on — most of
the frequency range was being thrown away as a single "not rare" bucket.
Switching analysis B to `lm(log10(max_af) ~ feature + PIP_qtl +
log10(cs_size_med))` uses the continuous MAF directly. On the real data this:

- raised the number of FDR-significant features from **5 to 7**,
- agreed in sign with the binary model's estimates on **87%** of features
  that were significant under either version,
- turned the promoter/enhancer contrast (PLS vs. dELS/pELS) from a
  non-significant trend under the binary cut into a decisive, FDR-significant
  result under the continuous model,
- was checked by permutation and remains valid despite non-normal residuals
  (heavy-tailed MAF), and
- lost exactly one feature relative to the binary version:
  `non_coding_transcript_variant`, whose association appears to be
  tail-specific (driven by the sub-1% MAF regime the binary cut isolated) and
  does not survive being pooled into the continuous scale.

**Direction convention.** `fit_feature_enrichment()`'s `family = "gaussian"`
branch fits the coefficient on log10(MAF); a negative coefficient means
carriers sit at *lower* MAF, i.e. rarer. That is transformed by `10^(-x)` (both
the point estimate and, after re-sorting, both CI bounds) so that, exactly
like analysis A's odds ratio, a value **above 1 means "carriers are rarer /
more enriched"** — both analyses read the same way on a shared log10 axis
with the reference line at 1. `.enrich_pairs()` tags the result
`estimate_type = "odds_ratio"` for A and `"maf_fold_change"` for B so the two
are never confused or plotted on an axis implying they're the same statistic;
`enrichment_forest_panel()` (`00_setup.R`) reads that column to switch the
supplement's B-facet x-axis label to "MAF fold change (carriers rarer →)"
while A keeps "Adjusted odds ratio".

`rare` remains a tested **feature** in analysis A. It is still excluded from
analysis B's feature list (`setdiff(feats, "rare")`) because it is
definitionally derived from the same `max_af` the continuous model now
regresses on directly — regressing a variable on a transform of itself is
meaningless, not because the binary label was removed from the pipeline
(the `rare` column, and its binary-model results, are no longer computed by
this analysis at all; there is no supplementary binary-B column retained).
Zero/negative `max_af` rows (which would give `log10(0) = -Inf`) are dropped
before fitting, with the dropped count reported via `message()`.

### `distTSS` / `gene_type` are gene-relative, not variant-level — split out of `.variant_feature_map()`

`.variant_feature_map()` originally folded `distTSS` and `gene_type` into its
one-row-per-variant map alongside genuinely variant-level columns (allele
frequencies, the cCRE flags, `consequence`, `phylop`), on the premise that all
of them are "properties of the variant alone". That premise is false for
these two: a variant's distance to a TSS and the type of gene it annotates
are properties of the variant **relative to a particular gene** — the same
variant can sit 14,508 bp from one gene's TSS and 13,376 bp from another's,
or annotate a protein-coding gene under one credible set and a lncRNA under
another. On the real staged data this hard-errored the whole variant-level
build: 180,211 variants carried "conflicting" values, all of them actually
correct per-gene rows collapsed into one map that could only hold one.

Diagnosed over a 100,000-variant sample of the real data: `gvs_*_af`, `PLS`
(and the rest of `FEATURE_CCRE_COLS`), `phylop`, `consequence` never carried
more than one distinct value per variant; `distTSS` did for 33,435 variants,
`gene_type` for 26,030.

Fix: `distTSS`/`gene_type` (`FEATURE_GENE_REL_COLS`) moved to a new
`.variant_gene_feature_map()`, one row per **(variant, gene)**, joined into
`coloc_variant_level()` on `c("SNP" = "variant", "Gene" = "gene")` — `Gene`
being the same bare-ENSG column `coloc_variant_level()` already derives from
the enloc `Signal`. This makes `tss_proximal`/`protein_coding` mean
"...for the gene this colocalization is about", which is what the analysis
intends; before the fix a shared variant could only ever get one arbitrary
gene's answer (when the guard didn't just abort the whole build).
`.variant_feature_map()`'s original guard is otherwise **unchanged** and
still fires on a genuine variant-level disagreement.

Verified against the real staged data: rebuilding `coloc_variant_level()`
now succeeds (0 variant-level conflicts, down from 180,211), with 268,437
variant-level rows and only 197 (0.07%) missing a gene-relative feature —
confirming the (variant, gene) join key is correct, not silently matching
nothing.

## Colocalization covariate analysis (dropped from the pipeline)

The original built a per-gene table named `coloc_summary` that was never
plotted or written. It joined, per gene:

- credible-set colocalization rate (`PropColoc`, `NumColoc`, `Total`) from the
  component collapse,
- median expression (`median_tpms.tsv.gz` → `mean_TPM`),
- GeneBayes s_het selection estimates (`post_mean`, plus quintile `shet_q` and
  the four `shet_bin` categories) from the Nature Genetics supplement
  `41588_2024_1820_MOESM4_ESM.xlsx`, sheet 2,
- protein detectability (`LOD_proportions_per_proteins.tsv`), joined via the
  gene → UniProt map,

filtered to genes with finite `mean_TPM` and `post_mean`.

The evident question: **does cross-assay colocalization rate track expression
level, selection pressure, or protein detectability?** The code is dropped;
this note is the record. To resume, rebuild it on top of
`coloc_cs_components()` and `susie_finemapping()$gene_map`.

Manifest keys `median_tpm`, `genebayes_shet` and `protein_lod` are kept for
provenance but have **no consumer in the current pipeline**. Note that
removing this analysis also removed both `left_join(genebayes, ...)` calls
that attached the s_het columns in the original — no retained panel reads
them.

Also dropped as unused: `cs_count`, `n_eQTLs`, `n_sQTLs`, `n_pQTLs`,
`rare_gene_list`, and the exploratory `glm(Coloc ~ gvs_max_af + PIP_qtl +
PIP_gwas_marginal)`.

## Open question — the eQTL-pQTL model formula

`panels/panelC.R` fits one model per assay pair (`FIT_FORMULA`). Two of them
use `abs()`:

```
eQTL-sQTL:  log1p(abs(posterior_mean.y)) ~ log1p(abs(posterior_mean.x)) * log_maf
sQTL-pQTL:  log1p(abs(posterior_mean.y)) ~ log1p(abs(posterior_mean.x)) * log_maf
eQTL-pQTL:  log1p(posterior_mean.y)      ~ log1p(posterior_mean.x)      * log_maf
```

The eQTL-pQTL formula omits it. This is preserved **verbatim** from the
original rather than silently unified, because it changes the estimate: SuSiE
posterior means are signed, and `log1p()` of a value below −1 is `NaN`, so any
eQTL-pQTL variant with a posterior mean below −1 is dropped from that fit but
retained in the other two.

If the asymmetry was unintended, add `abs()` to the `eQTL-pQTL` entry of
`FIT_FORMULA` in `panels/panelC.R`.

## Inputs still missing a `gs://` origin

`inputs.tsv` records these as local-only; they must be pre-staged, or given a
URI so the pipeline can self-stage on Terra:

| key | files |
|---|---|
| `susie_eqtl`, `susie_sqtl`, `susie_pqtl` | the three cleaned COMB SuSiE tables |
| `enloc_gene`, `enloc_cs`, `enloc_variant` | the nine FastENLOC outputs |
| `median_tpm`, `genebayes_shet` | deferred analysis above |

In the original these were hardcoded to
`/home/jupyter/multiomics_qtls/edit/FastEnlocInput/` or to bare
working-directory filenames.

## Fig3 reads the same SuSiE tables as Fig2

The `susie_eqtl` / `susie_sqtl` / `susie_pqtl` keys were repointed from

    finemapping/coloc/{eQTL,sQTL,pQTL}_susie.COMB.annotated.cleaned.tsv.gz

to the tables Fig2 already reads:

    COMB.{expression,splicing,proteomics}_SusieMerged.annotated.ancestry_skew.tsv.gz

**Verified against the real files.** All three carry every column Fig3 needs --
`cs_id`, `cs_index`, `posterior_sd`, `gene_name` (expression) and the seven
`AF_POP_COLS`. `susie_finemapping()` builds 844,442 eQTL / 254,361 sQTL / 60,636
pQTL rows, with 755 genes fine-mapped in all three assays.

The FastENLOC credible-set ids join back to these tables, which is the thing a
mismatched export would break silently:

| pair | edges | left matched | right matched |
|------|------:|-------------:|--------------:|
| eQTL-pQTL | 1,996 | 100% | 100% |
| eQTL-sQTL | 3,663 | 100% | 97.1% |
| sQTL-pQTL |   371 | 70.9% | 100% |

The sQTL side matches 97.1% in `eQTL_sQTL` but 70.9% in `sQTL_pQTL`. **This is
expected, not a defect.** Each assay pair's FastENLOC run is given its own
inputs by design, because the set of genes and credible sets shared between
eQTL/sQTL, eQTL/pQTL and sQTL/pQTL is different in each case. A pair's run
therefore covers only the credible sets relevant to that pair, so a lower match
rate against the full splicing table is what the design implies -- not evidence
that the wrong export was used.

The rates matter only as a check that the join CONNECTS at all: a repointed
SuSiE table whose `cs_id` format disagreed with the enloc output would show
near-zero matches on the affected side, which is the failure this table rules
out.

### The `group` filter is gone

The old pQTL read did `filter(group == "COMB")`. These are combined-ancestry
exports, one file per assay, so there is nothing to select on -- and the filter
would in fact have matched **nothing**: the `group` column holds
`"COMB.proteomics"`, not `"COMB"`. Left in place it would have handed the
pipeline a silently empty table.

The fixture previously carried a deliberately-wrong non-COMB pQTL row as a
canary for that filter being dropped; with the filter gone the row would be a
genuine duplicate, so it has been removed along with the assertions covering it.

### Schema validation

`.read_susie_table()` checks the schema up front and reports *every* missing
column at once, naming the offending key and file. A bare `select()` error names
only the first missing column, and these files come from an upstream pipeline
this repo does not control.

### Duplicate (variant, gene) keys are collapsed, not rejected

Splicing `molecular_trait_id`s are **intron-level**
(`chr7:150803456:150803620:clu_96688_+:ENSG00000002933.9`) and are collapsed to
a gene in `susie_finemapping()`. A variant sitting in credible sets for two
introns of the same gene therefore produces two rows for one (variant, gene)
key. In the shared tables:

| assay | rows | duplicated (variant, gene) | max copies |
|-------|-----:|---------------------------:|-----------:|
| sQTL  | 254,361 | 13,950 | 5 |
| eQTL  | 844,442 | 3 | 2 |
| pQTL  |  60,636 | 3 | 2 |

2,257 of 5,761 splicing genes have more than one intron cluster. The
`.cleaned` export this replaced was evidently pre-deduplicated, which is
presumably what "cleaned" meant; on the shared tables the many-to-one guard in
`coloc_model_input()` fired and stopped panels C and D from building.

`posteriors()` now keeps the row with the largest `|posterior_mean|` per
(variant, gene). Panels C/D regress cross-assay **effect sizes**, so the
strongest signal is the faithful single representative of what a variant does to
a gene; keeping both rows would double-count it and inflate the regression's N.

Two deliberate scoping decisions:

- The collapse happens at the **posterior join**, not in `susie_finemapping()`.
  The full tables still feed `cs_list`, `gene_map` and the credible-set graph
  behind panels A and B; dropping rows there would delete whole credible sets
  from them.
- The `relationship = "many-to-one"` guard is **retained**. It is no longer
  reachable from a duplicated input row, but it now asserts the collapse did its
  job. The tests that used to drive that error have been rewritten around the
  new contract.

**The 3 eQTL / 3 pQTL duplicates above are not all harmless.** Discovered
while verifying the `distTSS`/`gene_type` fix above against real data: for 5
of these 6 rows the duplicate copies agree on every variant-level column, so
`.variant_feature_map()`'s own `distinct()` silently collapses them for
free. For exactly one --
`chr3_23947795_C_CT` / `ENSG00000224074` in the eQTL table -- the two copies
disagree on `CTCF.bound` (`TRUE` vs `FALSE`) because they matched two
different ENCODE cCRE records (`ENCODE_ID_1`/`ENCODE_ID_2` differ between the
copies). Left alone this tripped `.variant_feature_map()`'s variant-level
guard on real data even after the gene-relative split, on a row that has
nothing to do with genes at all.

`.variant_feature_map()` now collapses duplicate `(variant, molecular_trait_id)`
rows **within each single SuSiE table**, first-row-wins, before binding the
three tables together -- the same "collapsed, not rejected" policy as
`posteriors()` above, applied to the same upstream duplicate-row class. The
collapse is scoped to within one table specifically so it can never absorb a
*genuine* disagreement between two *different* assay tables about the same
(variant, gene), which remains exactly what the guard exists to catch.

---

## Panel E — replaced: symmetric background + rare-variant enrichment

**This changes results.** Panel E previously rendered the pair-contrast scatter
of the colocalizing-vs-layer-specific feature enrichment (Analysis A of
`docs/superpowers/specs/2026-07-27-fig3-panelE-coloc-enrichment-design.md`).
It now shows, among **colocalizing eQTL-pQTL variants only**, which annotations
mark rarer and which mark commoner variants. The Analysis A machinery
(`coloc_feature_enrichment()`, `.enrich_pairs()`, `enrichment_forest_panel()`)
is untouched and still drives `Fig3/supp_coloc_enrichment.R`.

New code lives in **`Fig3/symmetric.R`**, sourced from `00_setup.R`.

### 1. The background was one-sided (the reason for the rewrite)

FastENLOC enumerates credible sets from ONE layer. In `Signal`, the half before
`(@)` is the QTL/annotation side; the half after it is the partner layer and is
populated only on colocalization. So every non-colocalizing row of
`.snp.out`/`.sig.out` is an annotation-side-specific credible set, and the
partner layer's unshared signals are **absent from the output entirely**.

For eQTL-pQTL the background was therefore 11,065 eQTL-specific credible sets
and **zero** pQTL-specific ones. `symmetric_universe()` reconstructs the missing
side from the SuSiE table at the genes each pair was tested on:

| pair | annotation side | partner side recovered | unrecoverable |
|---|---|---|---|
| eQTL-pQTL | 13,295 eQTL CS (2,230 coloc) | 9,768 pQTL CS (2,143 coloc / 7,625 specific) | 0 / 2,143 |
| sQTL-pQTL | 770 sQTL CS (390 coloc) | 2,338 pQTL CS (379 / 1,959) | 7 / 386 |
| eQTL-sQTL | 37,606 eQTL CS (4,515 coloc) | 26,800 sQTL CS (3,307 / 23,493) | **1,075 / 4,382** |

**Panel E uses eQTL-pQTL only.** It is the only pair whose partner side is
fully recoverable; eQTL-sQTL loses 24.5% of its colocalizing credible sets to
introns absent from the staged splicing table.

### 2. `10^coef`, not `exp(coef)`

The outcome is `log10(max_af)`, so the fold change is `10^coef`. The
exploratory code used `exp()`, which understated every effect (`stop_gained`
3.18x -> 7.62x, `dELS` 0.88 -> 0.75). `exp()` remains correct for Analysis A,
which is logistic. Pinned by a test.

### 3. Coding consequences redefined as "clean canonical"

`consequence` is a VARIANT-level VEP call on whichever transcript VEP picked.
`missense_variant`, `synonymous_variant`, `stop_gained` and `frameshift_variant`
now require the variant to be in a **canonical CDS of the eGene itself**
(GENCODE v48), with SpliceAI max delta <= 0.2 and no cCRE overlap. The
pre-cleaning column survives as `<feature>_raw`.

`splice_region_variant` is deliberately NOT cleaned: the SpliceAI criterion
discards exactly the variants that make the class meaningful (it cut n from 184
to 23 and collapsed the evidence), and leaving it raw makes it consistent with
`splice_donor`/`splice_acceptor`, which were never in the list.

This adds **`gene_gtf`** to `Fig3/inputs.tsv`.

### 4. One row per variant, PIP-weighted

~31% of colocalizing variants sit in BOTH layers' credible sets and were being
counted twice as independent observations. `collapse_to_variants()` collapses
to one row per variant weighted by the **geometric mean** of its PIPs, which is
high only when the variant is confidently fine-mapped in both layers.

`log10(cs_size)` was dropped from the model: `sum(PIP)` per credible set is
0.997 (median, IQR 0.967-1.000), so PIP weighting already normalises each signal
to one unit of evidence and cs_size was a second proxy for the same thing.
Dropping it changed no estimate beyond the third decimal and no significance
call.

**The PIP weighting is load-bearing and must be stated, not assumed.**
Unweighted, only `conserved` survives; `dELS` goes 0.80 -> 0.98 (null). The
weighted result is corroborated independently by the PIP > 0.5 subset (dELS
0.78). Weighting is a design decision about the estimand -- per signal rather
than per credible-set member -- not a tuning choice.

### 5. What the panel does NOT show

The associations are **descriptive of colocalizing variants**. They are not
colocalization-specific: `stop_gained` behaves identically among layer-specific
variants (ratio-of-fold-changes 1.31, p = 0.33), and `conserved` and `PLS` are
*weaker* among colocalizing variants than layer-specific ones. `splice_region`,
`3' UTR` and `dELS` do differ. "Rare colocalizing variants are enriched for X"
is supportable; "colocalization selects for X" is not.

### 6. Two silent parse traps, both now tested

Both returned an EMPTY join rather than an error:

- Both halves of `Signal` carry an `ENSG<digits>:` prefix the SuSiE `cs_id`
  lacks. Leaving it on matches nothing.
- Multi-CS signals split on `_` before a NEW gene-prefixed id. Splitting on a
  plain `_ENSG` shreds the UniProt-prefixed pQTL ids
  (`Q6PCB0_ENSG00000179403.12_L2`) into fragments.

A third, in the panel builder: a CI end outside the x scale is dropped by
`geom_linerange` with only a warning, so the limits are derived from the data.
And an annotation missing from `PANEL_E_LABELS` is now **named in a message**
rather than vanishing -- that is how `stop_gained`/`frameshift` disappeared
during development.

---

## Panel B — the tallest bar was being deleted, not clipped

**This changes the figure.** `panels/panelB.R:19-20` capped the y scale with
`scale_y_continuous(limits = c(0, 0.55))`. That is not a zoom: `limits =`
**removes** out-of-range values before drawing, so any bar above 55% is dropped
entirely rather than truncated.

The cap was inherited from the original notebook export (`04128bc`), where it
presumably matched the data of the day. Once Fig3 moved to Fig2's shared SuSiE
tables (`aed63c2`) the distribution shifted and the "1 colocalized signal" share
rose to **66.0% (2,298 of 3,483 genes)** — so the single largest category
vanished from the panel, leaving a bare axis tick at 1.

It was invisible for two compounding reasons: ggplot reports the loss only as a
"Removed 1 row containing missing values" warning among many, and panel B hides
its y axis (`axis.text.y = element_blank()`), so a truncated bar cannot be read
as truncated — it just isn't there.

Replaced with `expand = expansion(mult = c(0, 0.05))` and no hard limit, so
every bar is drawn at its true height. The percentages are unchanged; only the
rendering was wrong.

**The general rule, which bit three times in panel E as well:** use
`coord_cartesian(ylim = )` to zoom and `scale_*_continuous(limits = )` only when
dropping data is genuinely intended. `build_panelE_maf()` derives its x limits
from the data for exactly this reason.
