# Figure 4 — draft legend

> Draft legend generated from the panel code. Values marked `[TK: ...]` must be
> filled in from the rendered figure or the data — they are not recoverable from
> the code alone. Verify against the final rendered figure before submission.

## Main figure

**Figure 4. Colocalization of molecular quantitative trait loci with GWAS traits
and the contribution of rare variants to gene-level trait association.**

**(A)** Colocalization of molecular quantitative trait loci (QTLs) with
genome-wide association study (GWAS) traits, summarised by trait category. Each
row is one GWAS trait; traits are blocked by category (category label at the left
of each block, marked by a grey rule) and ordered within a block by decreasing
gene count, with categories themselves ordered by their total gene count.
Individual trait names are suppressed on the *y* axis. Left sub-panel
(*Colocalizing Genes*): the number of distinct genes with at least one
QTL–GWAS colocalization for that trait; bar shading is scaled within each
category from its lowest to its highest count and is redundant with bar length,
while hue encodes trait category. Right sub-panel (*CS coloc rate*): the
credible-set (CS) colocalization rate for the same trait, expressed as a
percentage — [TK: numerator and denominator of the "lenient union" credible-set
colocalization rate], taken from the `union` layer of the colocalization-rate
table, i.e. pooled across molecular-QTL layers [TK: which layers the union
pools]. The dashed grey vertical line marks the median CS colocalization rate
across the traits displayed. Only traits with more than five colocalizing genes
are shown, at most the ten highest-count traits per category, and the
neuropsychiatric and renal/urological categories are omitted; [TK: number of
traits and number of categories displayed]. **(B)** Distribution across genes of
the increase in transcriptome-wide association study (TWAS) model accuracy
contributed by rare variants. For each gene, the prediction-model *R*² was
normalised by that gene's estimated heritability (*H*²) [TK: what *H*² denotes in
`R2_TWAS_combined.tsv`], and the plotted quantity is the difference between the
normalised *R*² obtained with no allele-frequency (AF) filter (AF threshold 0)
and that obtained after restricting to variants with AF ≥ 0.01, expressed as a
percentage; positive values indicate accuracy gained by admitting rare variants
[TK: confirm the semantics of the `AF_threshold` column]. Genes were restricted
to those with *P* < 1 × 10⁻⁵ [TK: which test the `Pval` column reports]. The *y*
axis is the number of genes falling in each of 30 equal-width bins; the *x* axis
is truncated at −20% and +40%, so genes outside that window are not drawn
[TK: total number of genes, and number falling outside the display window].
**(C)** Quantile–quantile (Q–Q) plots of conditional-TWAS association *P* values
across nine blood traits (eosinophil, haematocrit, lymphocyte, mean corpuscular
haemoglobin (MCH), mean corpuscular haemoglobin concentration (MCHC), monocyte,
neutrophil, platelet count and red-cell distribution width (RDW)) [TK: whether
the cell-type traits are counts or percentages], shown separately for the
common-variant test (*Common*) and the rare-variant test conditional on the
common-variant signal (*Rare*). Genes that colocalize with the trait are compared
with a matched null set of non-colocalized genes, matched 1:1 by nearest-neighbour
propensity score (logistic-regression propensity on log(1 + mean transcripts per
million (TPM)), log₁₀ of the selection coefficient *s*<sub>het</sub>, log(1 +
number of rare variants) and |Vcor| [TK: definition of Vcor]), exact-matched on
trait, using a 0.2 standardised-deviation caliper and matching without
replacement; the random seed is fixed in the code so the matched set is
reproducible. Before matching, one record was retained per gene–trait pair (that
with the highest gene-level colocalization probability, breaking ties on the
smallest conditional rare-variant *P*), and putative artefacts — records whose
conditional rare-variant *Z* exceeded three times its unconditional value while
|Vcor| > 0.3 — were removed. Axes are the expected against the observed
−log₁₀(*P*), with the expectation taken as rank/(*n* + 1) within each test and
group; the dashed line is the *y* = *x* null expectation. Observed values are
capped for display at 40 (*Common*) and 12.5 (*Rare*). The colour key is
suppressed in the panel: colocalized genes are drawn in blue (#275D8E) and the
matched null set in grey [TK: confirm the rendered colour of the "Null" series;
see Open items] [TK: number of colocalized and matched genes plotted].
**(D)** Enrichment of rare- relative to common-variant conditional-TWAS signal
across four clinically motivated gene sets. For each gene set, a logistic
regression was fitted over colocalizing genes, with gene-set membership as the
outcome and the gene's maximum absolute conditional rare-variant *Z* score,
maximum absolute common-variant *Z* score, log₁₀ *s*<sub>het</sub> and log₂ TPM
as predictors. Points are the odds ratios (ORs) for the rare (coloured) and
common (grey) *Z*-score terms, each expressed per one-unit increase in |*Z*|;
horizontal bars are the corresponding 95% profile-likelihood confidence
intervals; the dashed vertical line marks OR = 1, and gene sets are ordered by
OR. The outcomes are: *Mendelian Gene (Any)*, genes with a Definitive or Strong
ClinGen gene–disease validity classification; *Mendelian Gene
(Immune/Hematological)*, the subset of those genes whose ClinGen disease name
matches a hand-curated haematological/immune keyword list (see `NOTES.md`);
*Druggable gene*, genes with any recorded drug–gene interaction; and *Approved
Drug Target*, genes with an interaction involving an approved drug. The analysis
universe is genes colocalizing with at least one trait, after the same artefact
filter as in (**C**) and after excluding genes without an *s*<sub>het</sub>
estimate [TK: number of genes entering each regression; number of genes in each
gene set].

Error bars are shown only in (**D**), where they are 95% profile-likelihood
confidence intervals on the odds ratio; no interval or error bar is drawn in
(**A**), (**B**) or (**C**), and no significance test is annotated on any panel.
Regression *P* values in (**D**) are computed but not displayed. In (**A**),
colour identifies trait category, abbreviated as: AutoImm, autoimmune; Inflam,
inflammatory biomarker; Allerg, allergic/respiratory; Derm, dermatological;
CardioM, cardiometabolic; Lipid, lipid/metabolic; Metabolic, metabolic;
ClinChem, clinical chemistry; Heme, haematological; Anthro, anthropometric;
Renal, renal; GI, gastrointestinal; ReproEnd, reproductive/endocrine; MSK,
musculoskeletal; Cancer, cancer; Neurodeg, neurodegenerative; Neurol,
neurological; Ophth, ophthalmological (which of these appear depends on the
data; NeuroPsy, neuropsychiatric, and Renal/Urol, renal/urological, are excluded
by the panel filter). Abbreviations: QTL, quantitative trait locus; GWAS,
genome-wide association study; CS, credible set; AF, allele frequency; TWAS,
transcriptome-wide association study; TPM, transcripts per million;
*s*<sub>het</sub>, heterozygous selection coefficient; OR, odds ratio; Q–Q,
quantile–quantile. Molecular QTLs and their colocalization with GWAS traits were
derived from *All of Us* Research Program data [TK: GWAS resource(s) whose
summary statistics were colocalized, and the source of the trait-category
assignment]; colocalization was called upstream of this figure [TK:
colocalization method and the posterior-probability threshold defining a
colocalization]. Gene sets in (**D**) draw on ClinGen gene–disease validity
evidence and a drug–gene interaction table [TK: resource versions and access
dates].

## Open items

- `[TK: numerator and denominator of the "lenient union" credible-set
  colocalization rate]` — the panel plots the `pct_lenient_union` column of
  `coloc_rate_by_trait.tsv` as-is (`RATE_COLUMNS` / `coloc_by_trait()` in
  `00_setup.R`). Nothing in this repository defines how "lenient" differs from
  "stringent", or what the rate is a percentage *of*; it is set upstream.
- `[TK: which layers the union pools]` — only `layer == "union"` rows are kept;
  the set of molecular layers (expression, protein, splicing) the union spans is
  not stated in the code or the manifest.
- `[TK: colocalization method and posterior-probability threshold]` — `any_coloc`
  arrives as a precomputed flag in `harmonized_coloc.gene.tsv.gz`; the code only
  coerces it to logical. A `GLCP` (gene-level colocalization probability) column
  is used as a sort key in `prepare_qq_match_source()`, but no method name or
  cutoff appears anywhere in the repository, so the method must not be asserted
  from the code.
- `[TK: GWAS resource(s) and trait-category assignment]` — `trait` and
  `trait_category` are columns of the input table; no GWAS source is named in
  `inputs.tsv`, `00_setup.R` or the panels.
- `[TK: number of traits and categories displayed in (A)]`, `[TK: total number of
  genes in (B)]`, `[TK: number falling outside (B)'s −20% to 40% window]`,
  `[TK: number of colocalized and matched genes in (C)]`, `[TK: number of genes
  entering each (D) regression and the size of each gene set]` — all counts
  depend on the controlled-access data, which is not present locally.
- `[TK: what H² denotes]` and `[TK: which test Pval reports]` in
  `R2_TWAS_combined.tsv` — `panelB.R` uses both columns (`adj_R2 = r.squared /
  H2`; `filter(Pval < 1e-5)`) without defining them; a `H2_upper` column is
  computed there and never used.
- `[TK: semantics of AF_threshold]` — the code pivots on the literal values `0`
  and `0.01` and subtracts, which implies "all variants" versus "AF ≥ 1%", but
  the column is not documented.
- `[TK: definition of Vcor]` — used in the artefact rule (`|Vcor| > 0.3`) and as
  a matching covariate, but never defined.
- `[TK: whether the nine blood traits are cell counts or percentages]` — the
  `BLOOD_TRAITS` vector gives only short names.
- `[TK: rendered colour of the "Null" series in (C)]` — this is a probable code
  defect, not just an unknown. `plot_qq()` maps colours by the literal level
  names `"Matched non-colocalized genes"` and `"Colocalized genes"`, but
  `panelC.R` renames the matched group to `"Null"` *before* calling `plot_qq()`.
  `"Null"` is therefore unmatched by the manual scale and will render in
  ggplot2's `na.value` grey rather than the intended `#2B2B2B`. Because the panel
  also sets `legend.position = "none"`, the legend text must define both colours.
- `[TK: ClinGen / drug-interaction resource versions]` — `inputs.tsv` records
  both as local-only with no `gs://` origin or version.
- **No ancestry stratification in Figure 4.** Unlike the other manuscript
  figures, no population-group code (AFR, AMR, EAS, EUR, SAS, MID) appears
  anywhere in `Fig4/`; the legend therefore makes no ancestry claim. Confirm this
  is intended before submission.
- **README-vs-code discrepancies.**
  - The top-level `README.md` (line 28) lists `Fig4/` as "Flat script — not yet
    cleaned". Fig4 is in fact the modular manifest-driven pipeline
    (`00_setup.R`, `panels/panel{A–D}.R`, `assemble.R`, `run_all.R`); the
    top-level README is stale.
  - `Fig4/README.md` and the `Fig4.ipynb` markdown describe (**A**) as the
    "credible-set coloc rate" without stating that `coloc_by_trait()` selects the
    **lenient** union rate (`rate_type = "lenient"`), even though
    `build_plot_data()`'s own default is `"stringent"`. The code wins: the figure
    shows the lenient rate.
  - The `Fig4.ipynb` Panel A markdown says the panel is "restricted to the top-10
    **categories** with more than five colocalizing genes". The code restricts to
    **traits** with `count > 5` and keeps the top ten **traits within each
    category**, and separately drops two named categories. The notebook wording
    is wrong.
  - Neither README nor NOTES mentions (**B**)'s `xlim(-20, 40)` truncation, which
    silently drops genes from the displayed histogram.
