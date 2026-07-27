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

## Panel E — reserved slot (undecided)

The bottom row of Figure 3 reserves a slot beside panel D. Its content has not
been decided. `panels/panelE.R` renders a blank so the layout and
`assemble.R` stay unconditional.

Mechanics for whoever fills it: the bottom row is `(D | E)` with
`plot_layout(widths = c(1.3, 1))`, so E gets roughly 43% of an 18" canvas —
about 7.7 × 6 inches. Replace the body of `panels/panelE.R` and nothing else
needs to change.

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
