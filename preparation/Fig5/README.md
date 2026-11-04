# Figure 5 data preparation

These scripts contain the current Figure 5 transformations. They do not contain
input data. Run them from the repository root with the shared preparation runner.
Set the local `paper`, `source`, `input`, and `susie` roots as described in the
preparation guide. All files are written below `QTL_PREP_ROOT`.

Run `fig5-select` before `fig5-main`. The first step reads the three original
SuSiE tables and the original Open Targets ClinGen reference. It keeps rare rows
and the rows needed for the selected benign comparators. The second step builds
the A–E inputs and the statistical tests. It compares the new results with the
original accepted summaries through `prep_reference()`. This check always reads
the original reference file, even when a prepared overlay is present.

The selection rules are preserved:

- ClinGen evidence must be Strong or Definitive. A gene can have more than one
  inheritance mode. Each variant counts once across assays within each mode.
- A/B use the original maximum population AF < 0.01. They retain all ClinVar
  classes and add no PIP or effect cutoff. A uses the maximum frequency across
  AFR, AMR, EAS, EUR, MID, and SAS. Ties use that order. B assigns EUR when its AF
  is at least 75% of the maximum across those six groups.
- ClinVar classes use pathogenic, uncertain, then benign assertion precedence.
  A nonempty classification counts as annotated in B.
- C first retains unannotated, VUS, or P/LP rows. It then uses PIP > 0.9 and
  posterior mean < -0.3. The final gene step checks rare high-PIP variants of
  either sign at genes from the third stage. It requires a distinct pathogenic
  reference and an unclassified or VUS candidate.
- D preserves the original benign-comparator rule: exclude missense records,
  then prefer pQTL, eQTL, and sQTL in that order. E preserves the explicit
  comparator targets and the two later exclusions in the original scripts.
  Both tables contain the accepted posterior effects and current REM labels.
- Fisher tests use distinct variants. Subset tests use disjoint retained and
  excluded rows. The AD/AR union counts a shared variant once. Holm correction
  uses the original two-test, ten-test, and thirteen-test families.

`fig5-s28` rebuilds the historical 574-variant eQTL set. It retains the original
uncertain-or-empty ClinVar rule. It takes the strongest decreasing qualifying
target per variant and uses a strict expression-reduction threshold. This is a
different set from current A/B.

`fig5-s30` is a partial recovery. The original script contains supplied aggregate
counts for three ranks and a fixed IFT80 percentile. The step evaluates only
that recorded table/formula expression and checks all four display rows. It does
not recover the missing genome-wide ranking calculation or its input records.
The supplied counts are not copied into this repository.

The S31 notebook now reads the ACVRL1 fine-mapping table, GENCODE gene/exon
rows, public cCRE annotations, and supplied PromoterAI predictions directly.
It computes the locus range, merged exon intervals, causal windows and frequency
labels, then constructs both plots. Its visible helpers are in
`notebooks/current/plots/Fig5_locus.R`.

The separate archived `fig5-s31` runner still writes an RDS overlay and retains
its accepted D renderer. The notebook does not use that runner or those objects.
Fine mapping, PromoterAI prediction, motif scanning and external annotation
generation remain upstream inputs.

Source paths and SHA-256 values are in `provenance.json`. The R test checks
variant-versus-pair counts, inheritance overlap, ClinVar precedence, AF ties,
the inclusive EUR75 boundary, strict rarity/PIP/effect boundaries, the final
either-sign gene step, disjoint Fisher rows, strict reduction thresholds,
the S31 inclusive locus PIP boundary, and the locus renderer contract:

```sh
Rscript tests/test_preparation_fig5.R
```

Local execution matched the accepted current A–E tables, the Fisher-test tables,
both S28 threshold tables, all four S30 display rows, and the S31 prediction
windows. Source inputs are hashed before and after each step by the shared
runner. No cloud job or model fitting is part of these steps.
