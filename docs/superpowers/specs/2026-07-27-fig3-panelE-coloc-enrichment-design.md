# Fig3 panel E — colocalization feature enrichment — design spec

*2026-07-27. Fills the reserved panel E slot documented in `Fig3/NOTES.md`.*

## 1. Scope

Panel E of Figure 3 is currently a blank placeholder. This spec fills it with a
feature-enrichment analysis of **colocalizing versus QTL-layer-specific
variants**, plus a supplementary figure carrying the full result set.

Two questions are answered, each a separate analysis over the same feature set:

- **Analysis A — what makes a variant shared.** Within an assay pair, how do
  colocalizing variants differ from non-colocalizing ones?
- **Analysis B — what makes a shared variant rare.** Among colocalizing
  variants only, how do rare ones differ from common ones?

Analysis A is the panel's subject. Analysis B was added because rare variants
are already visibly over-represented among colocalizing variants (below), which
makes the follow-up question worth asking rather than speculative.

## 2. Locked decisions

| Decision | Choice | Rationale |
|---|---|---|
| Unit | variant × assay pair | The grain `coloc_variant_level()` already returns. |
| Contrast (A) | `Coloc == TRUE` vs `FALSE`, **within pair** | Background matched by construction — both groups are fine-mapped QTL variants in the same two layers. |
| Contrast (B) | **continuous log10 MAF**, within `Coloc == TRUE` | Same question — among shared variants, what distinguishes the rare ones — but without dichotomising a continuous variable at 1%. See §5a. |
| Pairs | tested independently, never pooled | eQTL-pQTL and eQTL-sQTL are the substantive comparison; pooling would hide it. |
| Features | all binary | One statistic (an odds ratio) so both figures read coherently. |
| Statistic | binomial GLM, adjusted OR | Colocalization requires good fine-mapping in *both* layers; a raw Fisher test conflates biology with detectability. |
| MAF adjustment | **none** | Rare-vs-common is itself a tested feature in A, and the outcome in B. Stated as a limitation, not an oversight. |
| Main panel | pair-contrast scatter | Shows the eQTL-pQTL vs eQTL-sQTL distinction directly rather than across dodged points. |

## 3. Baseline counts

From `coloc_variant_level()` on the staged data. These drive the power
decisions below and should be re-checked if inputs change.

| pair | non-coloc | coloc | rare+coloc | common+coloc | rare % non-coloc → coloc |
|---|---:|---:|---:|---:|---|
| eQTL-pQTL | 61,042 | 1,725 | 130 | 1,592 | 5.9 → 7.5 |
| eQTL-sQTL | 200,494 | 3,232 | 228 | 2,998 | 6.7 → 7.1 |
| sQTL-pQTL | 1,659 | 285 | 18 | 267 | 2.4 → 6.3 |

Rare variants are modestly over-represented among colocalizing variants in all
three pairs. `sQTL-pQTL` is an order of magnitude smaller than the others.

**Analysis B is still attempted for all three pairs.** It is not excluded for
`sQTL-pQTL` a priori — with 18 rare colocalizing variants essentially every
feature will fall below the count floor, so the underpowered flag in §7 fires
on its own. The criterion decides; the spec does not hard-code the outcome.

## 4. Feature set

All features are binary, derived per variant.

- **cCRE flags (7)** — `PLS`, `pELS`, `dELS`, `CTCF.bound`, `CTCF.only`,
  `DNase.H3K4me3`, `FANTOM5`. Already logical in the SuSiE tables.
- **VEP consequence terms** — one-hot expanded from the comma-separated
  `consequence` column.
- **Binarized covariates** — `rare` (`max_af` < 0.01, matching
  `coloc_model_input()`), `conserved` (`phylop` > 2), `tss_proximal`
  (|`distTSS`| < 2000 bp), `protein_coding` (`gene_type == "protein_coding"`).

The two cut-points are declared once as named constants in `Fig3/00_setup.R`
(`PHYLOP_CONSERVED <- 2`, `TSS_PROXIMAL_BP <- 2000`) and stated in the panel
caption. They are conventional defaults, not tuned — if a result turns on the
exact cut-point that is a finding about fragility, and belongs in `NOTES.md`.

`.expand_vep_consequences()` (written for Fig1 this session) **moves to
`common/`**. It would otherwise have two callers with divergent copies. Fig1's
`00_setup.R` sources it from the new location; behaviour is unchanged.

## 5. Statistic

Per (analysis, pair, feature), one model:

```
A:  glm(coloc    ~ feature + PIP + log10(cs_size), family = binomial)   # adjusted OR
B:  lm (log_maf ~ feature + PIP + log10(cs_size))                      # Coloc == TRUE only
```

Analysis A reports an adjusted **odds ratio**. Analysis B reports the fitted
coefficient as a **fold-change in MAF** (`10^coef`): "carriers of this feature
have N-fold lower minor-allele frequency". The two statistics are NOT
comparable and are never plotted on a shared axis — panel E is analysis A
alone, and the supplement facets them separately.

Reported: adjusted odds ratio `exp(beta_feature)`, profile CI, p-value, and
`n_carriers`. FDR correction is applied **across features within each
(analysis, pair)** — not across pairs, which are separate questions.

The `PIP + log10(cs_size)` adjustment is the point of the design. In A it
separates real biology from "this variant was simply easier to fine-map in both
layers". In B it earns its place independently: rare variants are harder to
fine-map, so without it any feature correlated with fine-mapping precision
would masquerade as a rarity signal.

## 5a. Why analysis B is continuous, not binary

The binary `rare` (MAF < 1%) cut was measured against the continuous
alternative on the real staged data before this was adopted.

**Power.** Significant features rose from 5 to 7 (eQTL-pQTL 2→3, eQTL-sQTL
3→4), with 87% direction agreement across the 23 tests. More important than the
count: the promoter/enhancer contrast crossed from trend to decisive —
`PLS` −0.153 (padj 0.098 → 0.0005), `dELS` +0.141 (0.098 → 0.0003), `pELS`
+0.093 (0.70 → 0.026).

**The cut was badly placed.** Among colocalizing variants the median MAF is
~20% and only 6-7% fall below 1%, so the binary split isolated a thin tail and
discarded the rest of the frequency information.

**Non-normal residuals are not a practical problem here.** Residuals are
decisively non-normal (Shapiro-Wilk p ≈ 1e-30 to 1e-43, skew −0.70 to −1.05,
hard ceiling at log10(0.5)). A 999-permutation null on the `PLS` coefficient
returned p = 0.001 and p = 0.043 against model p = 1.6e-04 and 4.0e-02 — well
calibrated. At n = 1,722-3,226 the CLT carries the coefficient estimates, and
Wilcoxon agrees independently.

**One real loss, retained deliberately.** `non_coding_transcript_variant` is
significant under the binary test (OR 2.13, padj 0.016) and not under log MAF
(padj 0.36) — it is enriched specifically in the extreme tail without shifting
the bulk. The binary result is therefore kept as a supplementary column rather
than discarded: the two tests answer different questions ("enriched below 1%?"
vs "shifts mean frequency?").

## 6. Guardrails

Three failure modes were hit in Fig1 panel H during this session. All are
designed out here, because each one produced a figure that rendered cleanly
while being wrong.

1. **Count floor with logged drops.** A feature is reported only with **≥10
   carriers in the smaller group**. Below that it is dropped *and named in a
   message*. Fig1 panel H silently returned infinite odds ratios from zero
   cells, and those points then vanished from the plot with no warning.
2. **VEP-native spelling.** Feature names match VEP exactly
   (`3_prime_UTR_variant`, not `_utr_`). A casing mismatch silently discarded
   two annotations from panel H via an equality join.
3. **Coverage assertion.** The number of features tested is compared against
   the number expected and the difference reported. A shrinking feature set
   must not look like a clean run.

## 7. Underpowered-pair handling

`sQTL-pQTL` is structurally absent from the main panel — that panel is a
two-axis scatter of eQTL-pQTL against eQTL-sQTL. Its treatment in the
**supplement** is data-driven rather than a fixed exclusion:

> A pair is flagged **underpowered** for a given analysis if fewer than half
> its features clear the count floor. Flagged pairs render in muted grey with a
> caption note.

On the counts in §3 this will almost certainly flag `sQTL-pQTL` for Analysis B
and may well pass it for Analysis A. Deciding per-analysis by criterion is more
honest than dropping the pair wholesale, and it re-evaluates itself if the
inputs change.

## 8. Outputs

**Main — `Fig3/panels/panelE.R`.** Pair-contrast scatter: one point per
feature, x = adjusted log-OR in eQTL-pQTL, y = same in eQTL-sQTL, diagonal
reference. On-diagonal features are generic shared-signal markers; off-diagonal
ones distinguish protein- from splicing-colocalization. Off-diagonal points
labelled, significance encoded by fill. Analysis A only.

The slot is `(D | E)` at `plot_layout(widths = c(1.3, 1))` — roughly
7.7 × 6 inches. Only the body of `panels/panelE.R` changes; `assemble.R` is
untouched.

**Supplement — `Fig3/supp_coloc_enrichment.R`.** Full forest of every surviving
feature for both analyses across all three pairs, faceted by functional group,
adjusted OR with CIs. Standalone, not built by `run_all.R`, matching
`supp_coloc_qc.R`.

## 9. Data flow

- `.variant_af_map()` → **`.variant_feature_map()`**, carrying annotation
  columns rather than only allele frequencies. Its existing one-row-per-variant
  guarantee is retained; the many-to-one join assertion in
  `coloc_variant_level()` continues to depend on it.
- New memoized builder **`coloc_feature_enrichment()`**, returning the tidy
  frame `(analysis, pair, feature, odds_ratio, conf.low, conf.high, p.value,
  padj, n_carriers, flagged)`. Cached via `build_or_cache` like its neighbours.
- One plotting function serves both figures, since both consume that shape.

## 10. Testing

Pure seams, testable without the real inputs:

- `.expand_vep_consequences()` — already has behaviour worth pinning: multi-term
  rows set every term, `.` produces no column, existing columns are never
  clobbered.
- **`fit_feature_enrichment(df)`** — the per-feature model fit, given a frame.

Assertions:

| Test | Expectation |
|---|---|
| Synthetic enriched feature | recovers OR > 1, CI excludes 1 |
| Null feature | CI includes 1 |
| Feature below the count floor | dropped, and named in the message |
| Coverage assertion | reports the shortfall when features are dropped |
| Underpowered pair | flagged when <50% of features clear the floor |
| Tidy shape | both analyses return identical columns |

The Fig3 fixture extends with annotation columns (cCRE logicals, a
`consequence` string with multi-term rows, phyloP, distTSS). It is already
manifest-driven, so no path changes are needed.

## 11. Known limitations

- **No MAF adjustment** in Analysis A. Rare-vs-common is a tested feature, so
  adjusting for it would condition away the effect being measured. A reviewer
  may reasonably ask for the stratified alternative.
- **`max_af` NA behaviour is inherited** from `coloc_model_input()`: a variant
  missing any one of the seven per-ancestry frequencies gets `max_af = NA` and
  drops out — 152 variants on current inputs. This is the same choice flagged
  as arguable in `Fig3/NOTES.md`; it is not re-litigated here.
- **Analyses A and B share a feature set but not a background**, so their odds
  ratios are not directly comparable to one another. They are plotted as
  separate facets, never on shared axes.

## 12. Out of scope

- Interaction / effect-modification testing (`coloc ~ feature * rare`).
  Considered and rejected: ~130 rare colocalizing variants leave interaction
  terms underpowered.
- Continuous features (phyloP, distTSS, pLI as raw values). Binarized instead,
  to keep one statistic across both figures.
- Any change to panels A–D, or to `assemble.R`.
