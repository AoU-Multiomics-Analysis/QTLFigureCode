# Supplemental Figure 1 — notes and deferred work

## Result-affecting change: the TF double log (Panel C)

`SuppFig1.R` applied `log2(x+1)` to the TF TPM columns in `TF_genes`
(`SuppFig1.R:173`), then applied `log2(·+1)` **again** in the boxplot
(`:184`), the `lm` (`:179`) and the paired `t.test` (`:192`) — so Panel C plotted
and tested `log2(log2(TPM+1)+1)` while the axis read `log2(TPM + 1)`. The cleanup
returns raw TPMs from `tf_expression_table()` and logs exactly once, so Panel C's
values and the reported lm / t-test change.

## Dropped: the isoform histogram

`SuppFig1.R:291–299` read `isoform_results`, a variable that is never defined (its
source path is commented out at `:22`), so the block errors. It is not carried
over. To reinstate it, add an isoform-results table to `inputs.tsv` and a small
`extra/` script that reads it and plots isoforms-per-gene.

## GO_annotations.csv is staged but not yet consumed

`go_annotations` is declared in `inputs.tsv` and can be read with
`fread(input_path("go_annotations"))`, but nothing in the current pipeline uses it
(the original only `head()`-ed it). It is kept as a hook; wire it into an analysis
when needed.

## `compute_mean_TPM` no longer emits `SE`

The original emitted `SE = sd(value)/n`, which is `sd/n` (not the standard error
`sd/sqrt(n)`) and was never plotted. The column is dropped.

## `extra/gene_type_detection.R` counts over all AoU genes, not the AoU∩GTEx set

Its AoU means come from `compute_mean_TPM(gene_annotations())` applied directly
to the AoU GCT (`aou_tpms`), matching the original's `AoU_mean_TPMs_v45` — a
mean computed over AoU alone, before any join to GTEx. It deliberately does
**not** use `merged_tpms()`, which inner-joins AoU to GTEx and would narrow the
gene universe to genes present in both cohorts, undercounting detection.

## Cleanup: `plot_gene_comparison`'s `point_size` now actually does something

The original hardcoded `size = 3` in the `geom_point()` call and silently
ignored its own `point_size` parameter. The cleanup wires `point_size` into
`geom_point(size = point_size)`. Every call site in this pipeline passes the
default (3), so this is a behavioral no-op for the assembled figure — recorded
here in case a future caller passes a non-default value and expects it to take
effect.

## Layout: the lncRNA panel's extra margin was dropped

The original's bottom row was `TF_gene_comparison | lncRNA_comparison & theme(plot.margin = margin(l = 2, t = 5, unit = "cm"))`. Because `&` binds tighter than `|`, that margin applied only to the lncRNA panel — a notebook artifact rather than a deliberate choice. `assemble.R` uses `P$C | P$D` with no extra margin, so panel D (lncRNA) sits slightly tighter than in the original. This is a layout-spacing difference only; no plotted or tested value changes.

## Retained extras (not in the assembled figure)

`extra/gene_type_comparisons.R` (snRNA, miRNA), `extra/rrna_expression.R`,
`extra/hemoglobin_genes.R`, and `extra/gene_type_detection.R` reproduce the
per-gene-type comparisons the original computed but did not assemble. Each is
standalone and saves its own `SuppFig1_extra_*.pdf`.
