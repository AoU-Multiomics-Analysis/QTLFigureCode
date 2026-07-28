# Figure 2 — draft legend

> Draft legend generated from the panel code. Values marked `[TK: ...]` must be
> filled in from the rendered figure or the data — they are not recoverable from
> the code alone. Verify against the final rendered figure before submission.

## Main figure

**Figure 2. Allele-frequency spectrum, ancestry restriction and dosage-like
effect sizes of fine-mapped causal variants across expression, splicing and
protein molecular QTLs in *All of Us*.**

**(A)** Cumulative distribution of the allele frequency of statistically
fine-mapped causal variants, shown separately for expression QTLs (eQTLs),
splicing QTLs (sQTLs) and protein QTLs (pQTLs). Variants were retained if they
had a SuSiE posterior inclusion probability (PIP) greater than 0.9 and a defined,
non-zero allele frequency. The *x* axis is the variant's maximum allele frequency
across *All of Us* genetic-ancestry sub-populations, on a log10 scale; the *y*
axis is the cumulative percentage of all such causal variants for that assay
whose allele frequency is less than or equal to the value on the *x* axis, so the
curve accumulates from the rarest variants upward and reaches 100% at the most
common. The shaded region marks the rare-variant interval used throughout the
figure (allele frequency < 0.01); the fraction of causal variants falling in it
is the height of each curve where it crosses 0.01. *n* = [TK: number of
fine-mapped causal variants per assay]. **(B)** Distribution of the
population-enrichment odds ratio of every fine-mapped causal variant (PIP > 0.9),
pooled across all three assays and de-duplicated to one row per variant because
the ancestry-skew statistic is assigned at the variant level and is therefore
identical across assays. Variants are grouped on the *x* axis by their primary
population — the ancestry sub-population in which the variant attains its highest
allele frequency — comprising AFR (African), AMR (Admixed American), EAS (East
Asian), EUR (European), MID (Middle Eastern), SAS (South Asian) and ADMIX (other
or admixed ancestry). Within each population, rare (allele frequency < 0.01) and
common (allele frequency ≥ 0.01) variants are plotted as separate, side-by-side
distributions, indicated by the paired labels beneath the axis. Points are
individual variants, the shaded area is a violin density, and the horizontal bar
is the median of the distribution; no confidence intervals are shown. The *y* axis
is the ancestry-enrichment odds ratio on a log10 scale, and the dashed horizontal
line marks the odds ratio of 5 used, together with a false-discovery-rate (FDR,
Benjamini–Hochberg) adjusted *P* < 0.05, to call a variant ancestry-enriched. The
numbers above each population are the counts of ancestry-enriched variants
assigned to that population. *n* = [TK: total number of fine-mapped variants
plotted in B]. **(C)** Absolute allelic fold change (|log2 aFC|), the magnitude of
the effect of the alternate allele on gene expression, for fine-mapped expression
causal variants (PIP > 0.9) grouped by their most severe annotated consequence.
Each variant was assigned a single term by taking the highest-ranked entry of its
Variant Effect Predictor (VEP) consequence annotations and its candidate
*cis*-regulatory element class; the terms retained are stop gained, frameshift,
splice region (all splice acceptor, donor, donor 5th-base, donor-region,
polypyrimidine-tract and splice-region terms collapsed into one category),
promoter-like signature (PLS), proximal enhancer-like signature (pELS), distal
enhancer-like signature (dELS), 5′ untranslated region (UTR) and 3′ UTR. Effect
sizes were capped at |log2 aFC| = 3 for display. For each consequence the
distribution is split into rare (darker shade) and common (lighter shade)
variants; the slab is the density of |log2 aFC|, the point is the median and the
two nested horizontal bars are the 66% and 95% quantile intervals. Consequences
are ordered from top to bottom by decreasing median |log2 aFC|, and the dashed
vertical line marks a two-fold change (|log2 aFC| = 1). The value to the right of
each consequence is the rare-minus-common difference in mean |log2 aFC|, taken as
the coefficient on rarity from a per-consequence linear model of |log2 aFC| on
rare/common status; asterisks denote the FDR-adjusted significance of that
coefficient (\*\*\* *P* < 0.001, \*\* *P* < 0.01, \* *P* < 0.05). *n* = [TK: number of
expression causal variants and number of consequence categories in C].
**(D)** Expected expression of a heterozygous carrier when a gene's independent
causal variants are combined cumulatively. For each gene, variants that decrease
expression (log2 aFC < 0, top) and variants that increase expression
(log2 aFC > 0, bottom) were ordered separately by decreasing magnitude of effect
and added one at a time; after including *k* variants the cumulative allelic fold
change is 2 raised to the summed log2 aFC, and the expected heterozygous
expression is the mean of the reference haplotype (1) and that cumulative fold
change, expressed as a percentage of the expression of a homozygous-reference
individual. The *x* axis is the number of independent variants included (1, 2, 3
or 4 or more); each point is one gene at that cumulative step, the shaded area is
a violin density, and the black bar is the median. The dashed line marks 100%,
i.e. no net change. The shaded horizontal bands mark the dosage-like windows:
"deletion-like", an expected heterozygous expression of 45–60%, comparable to
what is expected from a heterozygous deletion, and "duplication-like", 145–400%,
comparable to a heterozygous duplication. The percentage above each group is the
proportion of genes at that cumulative step reaching the dosage-like threshold
(≤ 60% for decreasing variants, ≥ 150% for increasing variants). Observations
with an expected heterozygous expression of 400% or more were excluded from the
display. *n* = [TK: number of genes in D]. **(E)** Enrichment of expression
outliers among individuals carrying a burden of QTL alleles, as a function of the
predicted percent change in expression. The *x* axis, shared with **(F)**, is the
lower bound of the predicted percent-change bin; the *y* axis is the outlier
enrichment odds ratio on a log2 scale, and the dashed line marks an odds ratio of
1 (no enrichment). Vertical bars are [TK: confidence level of the ci_low/ci_high
interval carried in the input table] confidence intervals. Underliers
(individuals whose expression falls below the outlier threshold for a gene) are
shown in red and overliers (expression above the threshold) in blue.
**(F)** Median number of protein-coding genes per individual in each predicted
percent-change bin. Bars are the median across individuals and error bars span
the 25th to 75th percentile (interquartile range) of the per-individual counts;
these are not confidence intervals. Shaded boxes highlight the deletion-like
(Del-Like) and duplication-like (Dup-Like) percent-change bins, and the bracketed
annotations give the number of genes in each window ([TK: exact quantity the "5
Genes" and "108 Genes" annotations count — they are hard-coded literals in the
panel script]). **(G)** Number of validated dosage-like events across gene
categories. Each row is a gene category; the red point is the count of expression
underlier (deletion-like) events and the blue point the count of expression
overlier (duplication-like) events, joined by a grey segment, with the counts
printed at each end. Rows are ordered by the total number of events. Categories
shown are [TK: the gene-category labels on the y axis of G].

Throughout, causal variants are those fine-mapped by SuSiE with a posterior
inclusion probability (PIP) greater than 0.9 within their credible set, and rare
and common variants are separated at an allele frequency of 0.01 (panels A–C).
Where a distribution is shown (panels B–D), the central bar or point is the
median; only panel C draws quantile intervals (66% and 95%), only panel E draws
confidence intervals, and only panel F draws an interquartile range — no other
panel shows an error bar. Multiple-testing correction, where applied (panels B
and C), is the Benjamini–Hochberg false-discovery rate. Abbreviations: QTL,
quantitative trait locus; eQTL, expression QTL; sQTL, splicing QTL; pQTL, protein
QTL; PIP, posterior inclusion probability; MAF, minor allele frequency; aFC,
allelic fold change; VEP, Variant Effect Predictor; PLS, promoter-like signature;
pELS, proximal enhancer-like signature; dELS, distal enhancer-like signature;
UTR, untranslated region; FDR, false-discovery rate; OR, odds ratio; AFR,
African; AMR, Admixed American; EAS, East Asian; EUR, European; MID, Middle
Eastern; SAS, South Asian; ADMIX, other or admixed genetic ancestry. All panels
were computed from SuSiE fine-mapping of molecular QTLs in the *All of Us*
Research Program cohort ([TK: cohort size and release version]); the underlying
individual-level data are controlled-access and available through the *All of Us*
Researcher Workbench under its Data Use Agreement.

## Supplemental figure — ancestry composition of enriched variants across the allele-frequency spectrum

**Figure S[TK: supplemental figure number]. Ancestry-enriched fine-mapped
expression variants across the *All of Us* allele-frequency spectrum.**
Fine-mapped expression causal variants (SuSiE PIP > 0.9) were divided into 30
bins of equal size by their minor allele frequency in the full *All of Us*
cohort, and each bin is plotted at the geometric mean of its allele-frequency
range (*x* axis, log10 scale). Stacked areas show the proportion of variants in
each bin that are ancestry-enriched — defined, as in Figure 2B, as an
ancestry-enrichment odds ratio greater than 5 with an FDR-adjusted *P* < 0.05 —
coloured by the ancestry group in which the variant is enriched (AFR, AMR, EAS,
EUR, MID, SAS and ADMIX), so that the total height of the stack is the overall
proportion of variants in that bin that are ancestry-enriched. The black line,
plotted on the same *y* axis, is the proportion of variants in each bin observed
at non-zero allele frequency in all seven ancestry groups, smoothed with a
centred five-bin rolling mean and then with a locally estimated
scatterplot-smoothing (LOESS) fit; no confidence band is shown. The white
vertical line marks an allele frequency of 0.5. *n* = [TK: number of fine-mapped
expression variants in the supplemental figure].

## Open items

- `[TK: number of fine-mapped causal variants per assay]` (A) — the counts depend
  on the controlled-access SuSiE tables; the code only applies the filters.
- `[TK: total number of fine-mapped variants plotted in B]` — after
  de-duplication to one row per variant; not computable without the data.
- `[TK: number of expression causal variants and number of consequence categories
  in C]` — the retained consequence set is fixed by `remove_list` /
  `consequence_order`, but which of those terms actually have variants (and how
  many) depends on the data.
- `[TK: number of genes in D]` — number of distinct `molecular_trait_id` values
  with a non-missing `log2_aFC`.
- `[TK: confidence level of the ci_low/ci_high interval carried in the input
  table]` (E) — `panelE.R` plots `ci_low`/`ci_high` straight from
  `QTLBurdenOutlierEnrichment.tsv`. Neither the level (95%? 89%?) nor the test
  that produced the odds ratio is set anywhere in this repository; it comes from
  an upstream script. The name of the test must be supplied before submission,
  since the legend should not claim one.
- `[TK: exact quantity the "5 Genes" and "108 Genes" annotations count]` (F) —
  both strings are hard-coded literals in `Fig2/panels/panelF.R` (lines 30 and
  34) with no accompanying computation, so their denominator and meaning cannot
  be recovered from the code.
- `[TK: the gene-category labels on the y axis of G]` — the labels come from the
  `label` column of `QTLBurdenOutlierCount.rds`. The code only reveals that one
  category contains the string "DrugInteraction" (relabelled "Drug Target") and
  that labels contain parenthesised sub-terms that are converted to line breaks.
- `[TK: cohort size and release version]` — no sample size or *All of Us* data
  release appears anywhere in the Fig2 code or manifest.
- `[TK: supplemental figure number]` — `supp_ancestry_stacked.R` writes
  `FigS_AncestryStacked.pdf`; the manuscript numbering is not encoded anywhere.
- `[TK: number of fine-mapped expression variants in the supplemental figure]`.

Ambiguities and code-versus-README discrepancies:

- **Panel D applies no PIP filter.** `01_prep_data.R` line 121 filters only on
  `!is.na(log2_aFC)`, so panel D is computed over every row of the merged
  expression SuSiE table, not the PIP > 0.9 subset. The README's framing (lines
  4–7) says the whole figure summarises variants at "SuSiE, PIP > 0.9". The
  legend above therefore describes D's variants as a gene's "independent causal
  variants" without a PIP threshold; if the merged table is already restricted to
  credible-set members this is defensible, but the threshold claim should not be
  carried into D's sentence until confirmed.
- **Panel A's *x* axis is labelled "MAF" but plots `gvs_max_af`**, the *maximum*
  allele frequency across ancestry sub-populations, not the cohort-wide minor
  allele frequency. The same field defines `rare` in `add_derived()`
  (`00_setup.R` line 72). The supplemental figure, by contrast, computes a true
  MAF by folding `gvs_all_af`. Either the axis label or the legend wording needs
  to be reconciled; the legend above says "maximum allele frequency across
  sub-populations" for A–C and "minor allele frequency" only for the supplement.
- **The panel C effect-shift labels are hard-coded with a leading "+"**
  (`01_prep_data.R` line 110: `paste0("+", round(estimate, 2), stars)`). A
  negative rare-minus-common shift would render as e.g. "+-0.12". Check the
  rendered figure; this is a display bug, not a legend issue, but it affects
  whether the legend can describe the annotation as "the difference".
- **The definition of the ancestry-enrichment odds ratio and its *P* value is not
  in this repository.** `gvs_odds_ratio` / `gvs_p_value` arrive pre-computed in
  the `*.ancestry_skew.tsv.gz` inputs; only the FDR correction and the OR > 5,
  *P*adj < 0.05 calling rule are set here (`01_prep_data.R` lines 44–46). The
  legend describes what the threshold does, not how the odds ratio was derived —
  that sentence must be added from the methods.
- **Panels E and F are stacked with E's *x* axis text and title suppressed**
  (`panelE.R` lines 20–21), which reads as a shared axis. But E's *x* is the
  continuous `focal_lower_bound` and F's is the discrete, reordered `pct_bin`
  factor, so the two axes are not guaranteed to align. Confirm on the rendered
  figure before the legend claims a shared *x* axis.
- **The supplemental *y*-axis title describes only the stacked areas.** The black
  line is a different quantity (proportion shared across ancestry groups) drawn
  on the same axis under the title "Proportion of ancestry-enriched variants".
- **"Shared across all continents" includes the `oth`/ADMIX group**
  (`supp_ancestry_stacked.R` lines 30–35), which is not a continental group. The
  legend says "all seven ancestry groups" rather than "all continents".
- **The README's repository layout is stale.** It describes a standalone
  `Fig2_clean.R` and a nested `fig2/` subdirectory (README lines 40–64); neither
  exists — the scripts sit directly in `Fig2/`. The README's input table (lines
  137–144) also lists `input_data/…` paths, whereas `inputs.tsv` uses bare
  filenames for the three SuSiE tables and `Fig2Data/…` for the three burden
  inputs.
- **Fig2 has no `NOTES.md`**, unlike Fig1 and Fig3–5, so there is no record of
  result-altering changes to cross-check this legend against.
- Panel D's shaded bands (45–60% and 145–400%) do not exactly match the
  thresholds used for the printed percentages (≤ 60% and ≥ 150%);
  `00_setup.R` lines 171–172 versus lines 195–197. The discrepancy is small but
  the legend states both numbers separately rather than implying they are the
  same window.
- Panels D and E/G use different colours for the same decrease/increase contrast
  (D: #B65C73 / #2A9D8F; E and G: #D44D5C / #046E8F). If the manuscript intends
  one colour convention for deletion-like versus duplication-like, this should be
  unified before the legend states a figure-wide colour code.
