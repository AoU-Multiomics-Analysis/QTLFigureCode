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

### 2. Panel D now excludes variants with no allele frequency, instead of silently binning them as "Rare"

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

### 3. `refresh = TRUE` does not cascade through the derived-data cache

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
