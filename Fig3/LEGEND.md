# Figure 3 — draft legend

> Draft legend generated from the panel code. Values marked `[TK: ...]` must be
> filled in from the rendered figure or the data — they are not recoverable from
> the code alone. Verify against the final rendered figure before submission.

## Main figure

**Figure 3. Genetic signals are frequently shared across the expression, splicing
and protein layers, and the degree of sharing depends on allele frequency.**

**(A)** Gene-level colocalization. For every gene fine-mapped in all three assays
— that is, genes carrying at least one SuSiE credible set in the expression
quantitative trait locus (eQTL), splicing quantitative trait locus (sQTL) and
protein quantitative trait locus (pQTL) tables (n = [TK: number of genes
fine-mapped in all three assays]) — we counted how many of the three assay pairs
(eQTL–pQTL, eQTL–sQTL, sQTL–pQTL) showed a gene-level colocalization. Bars give
the percentage of those genes (y axis) falling into each count from 0 to 3
(x axis); percentages sum to 100% across the four bars. The bracket annotates the
percentage of genes colocalized in at least one layer ([TK: percent of genes with
a colocalization in ≥1 layer], computed at build time, not hard-coded).

**(B)** Signal-level colocalization. FastENLOC reports colocalization between
individual credible sets ("signals"); colocalized credible-set pairs whose two
sides map to the same gene were treated as undirected edges, and the connected
components of the resulting graph — pooled over all three assay pairs — define
signals that are shared across molecular layers. Bars give the percentage of
genes (y axis) carrying each number of distinct shared signals (x axis). The
denominator here is genes with at least one shared signal ([TK: number of genes
with ≥1 colocalized signal]), so no zero category is shown, and it is not
restricted to the three-assay set used in (A). Inset: for each assay
individually, the proportion of that assay's credible sets that participate in a
shared signal, among genes fine-mapped in all three assays; assays are ordered by
that proportion.

**(C)** Dependence of cross-layer effect-size concordance on allele frequency.
For each assay pair we fitted an ordinary least-squares model of the second
assay's SuSiE posterior mean effect size on the first assay's posterior mean and
its interaction with −log10(minor allele frequency, MAF), over variants called
colocalized at the variant level. Points show the estimated interaction
coefficient and horizontal lines the corresponding 95% confidence interval; the
dashed vertical line marks a null effect of zero. Effect sizes enter the model as
log1p-transformed magnitudes, and −log10(MAF) is standardized (mean 0, s.d. 1)
across all three assay pairs pooled, so coefficients are on a common scale and
are comparable between pairs. MAF is the maximum, across the seven *All of Us*
per-ancestry allele-frequency estimates (African, AFR; Admixed American, AMR;
East Asian, EAS; European, EUR; Middle Eastern, MID; South Asian, SAS; and Other,
OTH), of the frequency folded to the minor allele; variants lacking any of the
seven frequencies were excluded.

**(D)** Concordance of expression and protein effect sizes, stratified by allele
frequency. Each point is a variant called colocalized between the eQTL and pQTL
assays; axes show the log1p-transformed SuSiE posterior mean effect of that
variant on gene expression (x) and on protein abundance (y) for the same gene.
Facets partition variants by MAF as defined above into common (MAF ≥ 0.05),
low-frequency (0.01 ≤ MAF < 0.05) and rare (MAF < 0.01) bins. Red lines are
ordinary least-squares fits with shaded 95% confidence bands; the grey dashed
line is the identity line y = x. Inset labels report the fitted slope β and the
number of variants N per facet ([TK: per-facet slopes and sample sizes]). Axes
are clipped to ±the 97th percentile of the absolute transformed posterior means
pooled across facets ([TK: rendered axis limits]), so a small number of extreme
points fall outside the plotted range.

**(E)** [TK: panel E content — the slot is reserved and currently renders blank
(`Fig3/panels/panelE.R`); either fill it or drop the letter before submission].

Colocalization throughout was assessed with FastENLOC, run separately for each
assay pair, at three nested resolutions: the gene-level colocalization
probability (GLCP) in (A), the regional colocalization probability (RCP) of a
credible-set pair in (B), and the SNP-level colocalization probability (SCP) of
an individual variant in (C) and (D). No fixed posterior cutoff was applied.
Instead, rows were called colocalized by controlling the Bayesian false-discovery
rate (FDR) at 5%: posterior probabilities were ranked in decreasing order, the
local FDR (1 − posterior) was averaged over each prefix, and the threshold was
set at the posterior of the last rank whose running mean remained at or below
0.05. This was done independently for every assay pair at every level, so the
implied posterior cutoff differs between pairs and between panels ([TK: realized
posterior cutoffs per assay pair and level, reported to the console at build
time]). Error bars and shaded bands appear only in (C) and (D) and are 95%
confidence intervals from the corresponding linear model; (A), (B) and the inset
in (B) show point estimates with no interval. No correction for multiple testing
was applied to the regression coefficients in (C) or (D) beyond the FDR control
used to define the colocalized variant set. Genetic data are from the *All of Us*
Research Program; fine-mapping is combined-ancestry (COMB) SuSiE output, and
credible set denotes a SuSiE 95% credible set with per-variant posterior
inclusion probability (PIP).

## Supplemental figures

**Figure S[TK: number]. Quality control of cross-assay colocalization
posteriors.** **(A)** Kernel density estimates of the FastENLOC posterior
probabilities, coloured by assay pair, for the gene-level (GLCP), credible-set
(RCP) and variant-level (SCP) statistics; densities are drawn semi-transparently
so overlapping pairs remain visible. **(B)** Number of colocalizations recovered
as a function of an arbitrary posterior cutoff. For each assay pair and each of
the three posterior statistics, the count of records exceeding a threshold of
0.01, 0.1, 0.5 or 0.9 is plotted on a log10 scale (y axis) against the threshold
(x axis). Note that these panels deliberately bypass the Bayesian FDR
thresholding used in the main figure, in order to show the sensitivity of the
discovery count to the cutoff ([TK: counts at each threshold]).

**Figure S[TK: number]. Proportion of genes colocalized, by assay pair.** Bars
give, for each assay pair, the proportion of gene-level FastENLOC records called
colocalized at a Bayesian FDR of 5%. The denominator is every gene tested for
that pair, so — unlike Figure 3A — genes are not restricted to those fine-mapped
in all three assays ([TK: proportion per assay pair]).

## Open items

- `[TK: number of genes fine-mapped in all three assays]` — panel A's
  denominator. Computed at runtime by `susie_finemapping()` as the intersection
  of the three assays' `molecular_trait_id` sets; not printed by panel A.
  `Fig3/NOTES.md` states 755 for the currently manifested SuSiE tables, but that
  figure is documentation rather than code output and should be re-derived.
- `[TK: percent of genes with a colocalization in ≥1 layer]` — panel A's bracket
  annotation. Computed and emitted to the console by `Fig3/panels/panelA.R`;
  read it from the build log or the rendered panel.
- `[TK: number of genes with ≥1 colocalized signal]` — panel B's denominator
  (`sum(nn)`), never printed.
- Panel B inset values — the per-assay credible-set colocalization rates are not
  printed; read them from the rendered inset if the legend needs to quote one.
- `[TK: per-facet slopes and sample sizes]` and `[TK: rendered axis limits]` —
  panel D's β/N labels and the 97th-percentile `coord_cartesian()` limits are
  computed from the data at build time.
- `[TK: realized posterior cutoffs per assay pair and level]` — the FDR
  procedure derives a different cutoff for each of the nine (pair × level)
  combinations. `add_fdr_threshold_label()` reports each one via `message()`;
  they are only recoverable from a real run.
- `[TK: counts at each threshold]` and `[TK: proportion per assay pair]` — the
  supplemental figures' numeric content.
- `[TK: supplemental figure numbers]` — the scripts write
  `FigS_CrossQTLColocQC.pdf` and `FigS_CrossQTLPropGeneColoc.pdf`; nothing in the
  repo assigns them manuscript S-numbers.
- `[TK: panel E content]` — `panels/panelE.R` renders a blank plot, but
  `assemble.R` still calls `plot_annotation(tag_levels = "A")` over five slots,
  so a tag "E" is drawn over empty space. Either the panel is filled or the
  letter must be removed before submission; `NOTES.md` records the slot as
  undecided.
- **Panel A axis label vs. what is computed.** The x axis reads "Number of
  Colocalized Molecular Layers", but `panelA.R` sums `Coloc` over the three
  *assay pairs*, so the quantity is the number of colocalized pairs (0–3), not
  the number of layers. These are not the same: a gene colocalized in two pairs
  that share an assay involves three layers, not two. The legend above describes
  what the code computes; consider relabelling the axis.
- **README vs. code — allele frequencies and the dropped annotation input.** The
  README's [Data & inputs](README.md#data--inputs) table still lists the SuSiE
  files under the paths
  `finemapping/coloc/{e,s,p}QTL_susie.COMB.annotated.cleaned.tsv.gz`, whereas
  `inputs.tsv` now points all three keys at Fig2's
  `COMB.{expression,splicing,proteomics}_SusieMerged.annotated.ancestry_skew.tsv.gz`
  exports. The standalone AoU annotation table has been removed from the
  manifest entirely and allele frequencies are now taken from the SuSiE tables
  via `.variant_af_map()`; the header comment above `build_or_cache()` in
  `00_setup.R` still describes "the AoU annotation table (the largest input in
  the repo)" as a shared upstream. Consequence for the legend: variant-level
  allele-frequency coverage is no longer genome-wide, and colocalized variants
  absent from every SuSiE table are silently dropped from panels C and D
  (`[TK: number of variant-level rows with no allele frequency]`, reported as a
  `note:` line by `coloc_variant_level()`).
- **Panel C model asymmetry.** The eQTL–pQTL model in `FIT_FORMULA`
  (`panels/panelC.R`) omits the `abs()` applied to both posterior means in the
  other two pairs, so eQTL–pQTL variants with a posterior mean below −1 yield
  `NaN` under `log1p()` and are dropped from that fit alone. `NOTES.md` flags
  this as possibly unintended. If it is retained, the legend's description of
  (C) as modelling effect-size *magnitudes* is accurate for eQTL–sQTL and
  sQTL–pQTL but not for eQTL–pQTL; if it is fixed, no legend change is needed.
- **Panel B's two denominators.** The main bars are computed over all genes with
  a shared signal, while the inset restricts to genes fine-mapped in all three
  assays. This is a real difference in the code, not an oversight in the legend;
  confirm it is intended before submission.
- The README describes panel C as "Effect of −log(MAF) on cross-assay effect-size
  correlation". The fitted quantity is an interaction coefficient in a linear
  model of one assay's effect size on the other's, not a correlation
  coefficient; the axis label in the code ("cross-ome correlation") inherits the
  same loose wording.
