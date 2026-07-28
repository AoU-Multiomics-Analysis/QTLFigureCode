# Figure 5 — draft legend

> Draft legend generated from the panel code. Values marked `[TK: ...]` must be
> filled in from the rendered figure or the data — they are not recoverable from
> the code alone. Verify against the final rendered figure before submission.

## Main figure

**Figure 5. Fine-mapped rare expression-QTL variants at ClinGen disease genes,
and sequence-model dissection of two causal regulatory alleles.**

**(A)** Rare variants at ClinGen disease genes, partitioned by the genetic-ancestry
group in which each variant attains its highest allele frequency. Variants were
restricted to genes with a Definitive or Strong ClinGen gene–disease validity
classification and to rare variants, defined throughout as a maximum allele
frequency below 1% across all ancestry groups (minor allele frequency, MAF < 1%).
Top: the number of such variants per ancestry group, faceted by the ClinGen mode
of inheritance of the gene (AD, autosomal dominant; AR, autosomal recessive;
genes curated only as semidominant, SD, are excluded from the facets). Each bar
is shaded by whether the variant carries a ClinVar classification (darker,
"ClinVar annotated") or not (lighter, "novel variant"). Ancestry groups are AFR
(African/African-American), AMR (admixed American), EAS (East Asian), EUR
(European), MID (Middle Eastern), SAS (South Asian) and ADMIX
[TK: preferred manuscript label for the group the code recodes from `oth` to
`admix`]. Bottom: for the same variants, the proportion of variants in each
inheritance class falling into each ClinVar category — pathogenic/likely
pathogenic (P/LP), variant of uncertain significance (VUS), benign/likely benign
(B/LB), or unannotated (no ClinVar record). Categories are assigned by matching
the ClinVar classification string in the order pathogenic → uncertain
significance → benign, so any classification containing "pathogenic" (including
conflicting interpretations) is counted as P/LP. Denominator of each proportion
bar is all rare ClinGen variants of that inheritance class shown in the panel
above. [TK: total number of rare fine-mapped ClinGen variants, and their
breakdown by ancestry, inheritance class and ClinVar category.] [TK: confirm
whether a posterior-inclusion-probability threshold was intended for this panel —
`build_panelA` applies none, unlike (B) and (C), which require PIP > 0.9.]
**(B)** Magnitude of the fine-mapped expression-QTL (eQTL) effect for rare
variants at disease and control genes. Variants were retained if they had a
posterior inclusion probability (PIP) greater than 0.9 in SuSiE (sum of single
effects) fine-mapping, were rare (MAF < 1%), and fell at either a ClinGen
Definitive/Strong gene or a control gene, the latter defined as a non-ClinGen
gene with a GeneBayes posterior-mean heterozygous selection coefficient (s_het)
below 0.01. The y axis is the mean of the absolute SuSiE posterior-mean effect
size across the variants in each group; it is not the allelic fold change.
Groups on the x axis are the aggregate of all retained variants, and the same
variants split by the direction of their allelic fold change (aFC) into
decreasing (log₂ aFC < 0) and increasing (log₂ aFC ≥ 0) variants; colour
distinguishes control genes, autosomal-dominant and autosomal-recessive ClinGen
genes. Error bars are 95% confidence intervals computed as the mean ± 1.96 ×
(standard deviation / √n) of the absolute posterior-mean effect within each
group; no formal statistical test is performed. Genes curated by ClinGen under
both AD and AR contribute their variants to both inheritance classes.
[TK: number of variants and mean effect ± 95% CI in each of the nine
category × direction cells.]
**(C)** Effect sizes of the strongest expression-decreasing rare variants at
ClinGen disease genes. For each gene, the single rare, confidently fine-mapped
(PIP > 0.9) variant with the most negative log₂ aFC was selected, and the seven
genes with the largest decreases are shown per inheritance class; HEXA and MLC1
are retained irrespective of rank. The x axis is the log₂ allelic fold change of
that variant (the log₂ ratio of expression from the alternative versus the
reference haplotype), restricted by construction to decreasing variants
(log₂ aFC < 0) and plotted over the range −7 to 0; the y axis lists gene symbols,
faceted by mode of inheritance. [TK: identity of the genes plotted and their
log₂ aFC values.] [TK: state whether the optional homozygote-free variant filter
was applied in the rendered figure — the code applies it only when the
controlled-access `NoHomozygoteVariants.txt` list is present, and ranks without
it otherwise.]
**(D)** Statistical fine-mapping of the *ACVRL1* expression QTL. Each point is a
variant tested in the *ACVRL1* cis window; the x axis is genomic position on
chromosome 12 (hg38, Mb) and the y axis the SuSiE posterior inclusion
probability. Point fill denotes the credible set to which a variant is assigned;
variants with PIP ≥ 0.9 are drawn as large diamonds and all other variants as
small circles, and rare (MAF < 1%) high-PIP variants carry an additional open
diamond outline. The dashed horizontal line marks the PIP = 0.9 threshold used
throughout the figure to call a variant causal. Vertical shaded bands mark ENCODE
SCREEN candidate cis-regulatory elements (cCREs) that overlap a high-PIP variant,
coloured by class (promoter-like signature, PLS; enhancer-like signature, ELS).
The boxed region is expanded in the inset, which spans ±500 bp around the
promoter variant at chr12:51,912,316; the track beneath the panel shows GENCODE
gene models over the same x range, with exons as filled blocks and transcription
start sites marked by arrows indicating strand. [TK: number of credible sets
resolved at *ACVRL1*, and the PIP, MAF and log₂ aFC of the lead promoter
variant.]
**(E)** Base-resolution PromoterAI predictions around the fine-mapped *ACVRL1*
promoter variants. Each facet is a ±20 bp window centred on one variant with
PIP > 0.9, annotated above with its genomic coordinate and alleles, its MAF and
the ancestry group in which that frequency is attained, and its log₂ aFC; the
diamond above each facet is coloured by the credible set of that variant. Points
give the predicted PromoterAI effect of every possible single-nucleotide
substitution at each position of the *ACVRL1* transcript ENST00000419526.6 within
the window, plotted against position relative to the fine-mapped variant. Points
with an absolute predicted effect of at least 0.5 are coloured by alternative
allele and all weaker substitutions are grey; the dotted horizontal lines mark
that ±0.5 threshold and the dashed line marks zero. The fine-mapped variant
itself is drawn as a large outlined circle. The y axis uses a signed-square
transform so that small base-level effects remain visible alongside the causal
substitution. The shaded box denotes an SP1 GC-box motif, annotated manually in
the code at −3.5 to +5.5 bp relative to chr12:51,912,316 C>T rather than derived
from a motif scan. [TK: number of causal *ACVRL1* promoter variants faceted, and
the PromoterAI score of the focal C>T substitution.]
**(F)** Statistical fine-mapping of the *SPR* expression QTL, plotted exactly as
in (D) over chr2:72,881,700–72,893,400 (hg38); the inset expands ±200 bp around
chr2:72,887,978. [TK: number of credible sets resolved at *SPR*, and the PIP, MAF
and log₂ aFC of the lead regulatory variant.]
**(G)** Borzoi sequence-model predictions for the *SPR* regulatory variant at
chr2:72,887,998 (header, with its MAF and log₂ aFC; the diamond is coloured by
its credible set). Tracks show predicted signal for the reference (filled area
and coloured line) and alternative (red line) alleles as a function of position
relative to the variant, in three progressively narrower windows: repressive and
broad marks (H3K36me3, EZH2, H3K27me3, H3K9me3) over ±3 kb, activating histone
marks and transcription (EP300, H3K27ac, H3K4me1, H3K4me3) over ±1.5 kb, and
CTCF, cohesin, chromatin accessibility (DNase) and CAGE over ±500 bp. Only marks
whose predicted signal changes by at least 10% between alleles are displayed, and
the predicted RNA track is omitted; the percentage beside each track label is the
change in mean predicted signal within ±96 bp of that mark's reference-track
peak, itself located within ±750 bp of the variant, expressed as
100 × (alternative − reference) / reference. The bracketed number at the top left
of each track is its data range. The orange band marks the ENCODE SCREEN element
containing the variant (EH38E2009446, proximal enhancer-like and CTCF-bound,
spanning −163 to +70 bp relative to the variant) and the vertical grey line marks
the variant itself. The bottom sequence logo shows base-resolution Borzoi
contribution scores for the predicted CTCF track from −8 to +12 bp around the
variant, with letter height proportional to the attribution of that base and the
positions altered by the variant highlighted in red. [TK: identity of the marks
retained in the rendered figure and their predicted percentage changes; the
Borzoi model version and the cell type or tissue context of the predicted
tracks.]

Throughout, rare denotes a maximum allele frequency below 1% across ancestry
groups, and a variant is called causal when its SuSiE posterior inclusion
probability exceeds 0.9. Error bars appear only in (B), where they are 95%
confidence intervals from the normal approximation (mean ± 1.96 × SEM); no
statistical test or multiple-testing correction is computed anywhere in the
figure. Control, AD and AR gene classes share a common colour scheme in (B) and
(C); ancestry groups share a common colour scheme in (A); and credible sets are
coloured by a categorical palette applied independently within each of (D)–(G),
so credible-set colours are not comparable between loci. Abbreviations: eQTL,
expression quantitative trait locus; QTL, quantitative trait locus; PIP,
posterior inclusion probability; aFC, allelic fold change; MAF, minor allele
frequency; s_het, heterozygous selection coefficient; AD, autosomal dominant;
AR, autosomal recessive; SD, semidominant; P/LP, pathogenic/likely pathogenic;
VUS, variant of uncertain significance; B/LB, benign/likely benign; cCRE,
candidate cis-regulatory element; PLS, promoter-like signature; ELS,
enhancer-like signature; TSS, transcription start site. Expression QTLs were
mapped and fine-mapped with SuSiE in the *All of Us* Research Program cohort
[TK: number of participants and the per-ancestry-group sample sizes of the eQTL
discovery set — not present anywhere in the Figure 5 code]; gene constraint
estimates are GeneBayes s_het posterior means, gene–disease validity assertions
are from ClinGen, clinical classifications from ClinVar, candidate
cis-regulatory elements from ENCODE SCREEN, gene models from GENCODE, promoter
variant effects from PromoterAI, and allele-specific regulatory predictions and
contribution scores from Borzoi.

## Supplemental figures

Three supplemental scripts live under `Fig5/extra/supplements/`. None is part of
the assembled Figure 5, and none is run by `Fig5/run_all.R`.

**Supplementary Figure (ADA allelic series;
`Fig5/extra/supplements/plot_ADA_supp.R` → `FigureS_ADA.pdf`). Rare-variant
allelic series at the *ADA* lymphoid enhancer.** **(A)** Fine-mapping of the
*ADA* expression QTL over chr20:44,643,000–44,651,000 (hg38); the x axis is
genomic position (Mb) and the y axis the SuSiE posterior inclusion probability,
with the dashed line at PIP = 0.9. The shaded band marks the enhancer at
approximately chr20:44,646,150–44,646,600. Five focal variants are drawn as
diamonds coloured by the direction of their effect on *ADA* expression and
labelled with their alleles, the ancestry group in which they are frequent, and
their log₂ aFC: chr20:44,643,792 C>T (common, −0.32), 44,646,218 C>T (SAS,
+0.96), 44,646,238 C>T (EAS, −1.00), 44,646,509 T>C (EUR, +0.65) and
44,650,274 C>A (AFR, −0.68); all other variants are grey. These five effect
sizes are hard-coded constants in the script rather than recomputed.
**(B)** Predicted immune-cell DNase accessibility from Borzoi for the reference
and alternative alleles of the paired ETS-site variants, over ±800 bp around each
variant; the shaded band marks ±60 bp around the variant. Titles quote
accessibility changes of +34% for 44,646,218 C>T and −50% for 44,646,238 C>T,
also hard-coded in the script. **(C)** Base-resolution contribution-score logos
over −7 to +9 bp around the same two variants, showing an ETS core created by the
alternative allele at 44,646,218 and disrupted at 44,646,238; the shaded box
marks the ETS core and the outlined column the variant position. **(D)** Heatmap
of the predicted change in immune transcription-factor binding (alternative
versus reference, as a percentage, clipped at ±60%) for the fourteen
most-affected factors across all five variants of the allelic series.
[TK: the identity of the transcription factors shown and their percentage
changes; the cell type or tissue context of the predicted immune tracks.]

**Supplementary Figure (ADA, journal-format variant;
`Fig5/extra/supplements/ada_scifig/` → `FigureADA.pdf`, 183 × 170 mm).** A
publication-formatted re-cut of the same *ADA* analysis with lower-case panel
labels. **(a)** Lollipop plot of the measured eQTL effect (log₂ aFC) of the five
allelic-series variants against genomic position on chromosome 20 (kb), coloured
by ancestry group, with the enhancer shaded. **(b)** Measured eQTL effect against
the predicted accessibility-attribution difference (alternative − reference) for
the same variants, with a dashed least-squares line fitted to the four
concordant variants; Pearson correlation coefficients are annotated for the four
concordant variants and for all five, and the Borzoi-non-concordant variant is
drawn as an open symbol. [TK: the two correlation coefficients — they are
computed from the data at run time.] **(c)** Reference and alternative
contribution-score logos over −6 to +8 bp for the two ETS-site variants, showing
core creation and disruption. **(d)** Predicted accessibility-attribution
difference per variant, ranked, with the measured log₂ aFC annotated on each bar.

**Supplementary Figures (per-gene Borzoi browsers;
`Fig5/extra/supplements/supp_browser.R` → `supp_<gene>.pdf`).** Allelic Borzoi
browser tracks, plotted as in Figure 5G, for the fine-mapped candidate variants
at *FAT1*, *SLC1A2*, *SLC6A19*, *SCN11A* and *BSND*. For each gene, reference
(filled area) and alternative (red line) predicted signal is shown for every mark
whose predicted signal changes by at least 8% between alleles; the shaded band is
a generic ±100 bp window around the variant rather than an annotated regulatory
element, and the title gives the gene, the variant coordinate and its alleles.
[TK: which of the five genes produced output, and the variant coordinates and
percentage changes for each.]

## Open items

Unresolvable from the code (require the rendered figure or the controlled-access
data):

- **(A)** Total number of rare fine-mapped ClinGen variants and their breakdown by
  ancestry group, inheritance class and ClinVar category.
- **(A)** Preferred manuscript label for the ancestry group the code recodes from
  `oth` to `admix`.
- **(B)** Per-cell sample sizes, mean absolute effects and 95% confidence
  intervals for the nine category × direction combinations.
- **(C)** Gene symbols shown and their log₂ aFC values.
- **(C)** Whether the optional homozygote-free variant list
  (`NoHomozygoteVariants.txt`, controlled access, `load_nohomozygotes()`) was
  present when the published figure was rendered. Its presence changes which
  variants are eligible, so the panel's content is not determined by the code
  alone.
- **(D)** Number of credible sets at *ACVRL1* and the PIP, MAF and log₂ aFC of the
  lead promoter variant.
- **(E)** Number of causal promoter variants faceted and the PromoterAI score of
  the focal chr12:51,912,316 C>T substitution.
- **(F)** Number of credible sets at *SPR* and the PIP, MAF and log₂ aFC of the
  lead regulatory variant.
- **(G)** Which marks survive the ≥10% filter in the rendered figure and their
  predicted percentage changes; also the Borzoi model version and the cell
  type/tissue context of the predicted tracks (the bigWigs are named by mark
  only, e.g. `CTCF.ref.bw`).
- **Cohort description.** Sample size and per-ancestry-group counts of the eQTL
  discovery set appear nowhere in the Figure 5 code; take them from the Methods.
- **Supplements.** Transcription factors and percentage changes in the ADA
  heatmap; the two correlation coefficients in `ada_scifig` panel (b); which
  per-gene browsers were produced by `supp_browser.R`.

Ambiguities and documentation that contradicts the code:

- **Top-level `README.md` (line 29) lists `Fig5/` as "Empty (planned)".** It is
  out of date: `Fig5/` is a complete seven-panel manifest-driven pipeline.
- **Panel A applies no PIP threshold.** `Fig5/README.md` and the design spec both
  describe panel A as showing "fine-mapped rare ClinGen variants", and the y-axis
  label in `build_panelA` reads "Fine-mapped Rare Variants", but the only filters
  in the code are ClinGen membership and MAF < 1%. Panels B and C both require
  PIP > 0.9. Either the SuSiE input table is already restricted to credible-set
  members (in which case the label is correct) or panel A counts all tested rare
  variants at those genes. This must be resolved before the legend is final; the
  computed-but-unused `causal = !is.na(log2_aFC)` flag in `build_panelA` suggests
  a filter was contemplated and dropped.
- **Genes curated only as semidominant are labelled "Control" in (B).**
  `variant_effect_by_inheritance` drops SD rows from `gene_assignments` before the
  join, so a ClinGen Definitive/Strong gene curated *only* as SD receives no
  inheritance assignment and is relabelled "Control", even though it is a disease
  gene and is not required to have low s_het. Confirm this is intended.
- **AD+AR double-counting is deliberate.** A gene curated by ClinGen under both AD
  and AR contributes its variants to both inheritance classes in (B) and (C).
  This is recorded in `Fig5/NOTES.md` and the design spec as a preserved
  behaviour matching the published figure, but it should be stated in the Methods
  if the counts are quoted.
- **Panel B plots the SuSiE posterior mean, not aFC.** The axis label reads
  "|Average QTL effect|", which is easy to misread as an allelic fold change;
  `Fig5/NOTES.md` flags the same point.
- **Design spec versus `assemble.R` (cosmetic).** The spec (§2, §6) specifies
  explicit `labs(tag=)` calls in `assemble.R`; the implementation instead uses
  `plot_annotation(tag_levels = list(c("A","B","C","D","F","E","G")))`. The
  rendered tags are A–G either way, and `Fig5/NOTES.md` describes the
  implementation correctly.
- **Unused manifest key.** `inputs.tsv` declares `promoterai`
  (`Fig5/inputs/promoterai_tss500.tsv.gz`), but no script reads it; panel E reads
  the committed `acvrl1_promoterai` table under `Fig5/reference/`. It does not
  affect the figure, but the manifest suggests a genome-wide PromoterAI input
  that is not used.
- **The SP1 motif in (E) is hand-annotated.** The motif box and its "SP1" label
  are supplied as a literal tibble in `build_panelE`, not derived from a motif
  scan; the legend should not imply that a motif-matching analysis was run.
