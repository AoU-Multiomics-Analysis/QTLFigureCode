# Figure 1 — draft legend

> Draft legend generated from the panel code. Values marked `[TK: ...]` must be
> filled in from the rendered figure or the data — they are not recoverable from
> the code alone. Verify against the final rendered figure before submission.

## Main figure

**Figure 1. Multi-omic cis-quantitative trait locus (cis-QTL) discovery across
genetic-ancestry groups in the *All of Us* Research Program.**

**(A)** Schematic of the study design [TK: describe the externally supplied
schematic — the panel is a blank slot in the code (`panels/panelA.R`) and its
content is not recoverable from this repository]. **(B)** Number of participants
in each genetic-ancestry group among the individuals included in the combined
multi-ancestry expression analysis set. Participants were assigned to groups by
the *All of Us* ancestry predictions and restricted to the research identifiers
in the combined (COMB) expression QTL (eQTL) sample list; the "other" prediction
is shown as ADMIX (genetically admixed / unassigned). Groups are AFR (African),
AMR (admixed American), EAS (East Asian), EUR (European), MID (Middle Eastern),
SAS (South Asian) and ADMIX; bars are ordered by decreasing participant count
[TK: n participants per ancestry group and the total]. **(C)** Proportion of
tested molecular phenotypes with a significant cis-QTL as a function of sample
size, for expression (left), and for splicing and protein abundance (right).
For each ancestry group and assay, cis-QTL mapping was performed with tensorQTL
permutations and the beta-approximated permutation *P* values were corrected
within that group–assay stratum by the Benjamini–Hochberg procedure; the *y*
axis is the fraction of all tested phenotypes in that stratum reaching a false
discovery rate (FDR) below 0.05. The *x* axis is the number of participants in
thousands on a log10 scale. Circles denote single-ancestry analyses and the
diamond denotes the combined multi-ancestry analysis (labelled AoU), whose
sample size is the sum of the per-group counts in (B); the parenthetical value
in the AoU label is the number of phenotypes in that analysis [TK: confirm what
the number printed beside the AoU point should represent — see Open items].
Points are coloured by ancestry group as in (B) [TK: number of phenotypes tested
per assay; sample sizes for the protein and splicing assays]. **(D)** Number of
genes with a cis-QTL in the combined multi-ancestry analysis that were not
already reported by external reference studies, stratified by gene biotype and
by assay. Genes were called using the tensorQTL permutation Storey *q* value
(*q* < 0.05) and were classed as novel if absent from the significant gene sets
of the corresponding references: expression genes (eGenes) versus GTEx v10 whole
blood (*q* < 0.05), INTERVAL and TOPMed whole blood (*q* < 0.05); splicing genes
(sGenes) versus GTEx v10 whole blood sGenes and INTERVAL cis-sQTLs; and protein
genes (pGenes) versus the cis-associations of the UK Biobank Pharma Proteomics
Project (UKBB-PPP). Biotypes were taken from GENCODE v48 and collapsed to
protein-coding, long non-coding RNA (lncRNA) and other non-coding RNA;
pseudogene, immunoglobulin, T-cell-receptor and TEC biotypes were excluded, and
lncRNA was excluded from the pGene comparison because the proteomic assay cannot
measure it. Bars are drawn on a log10 scale (counts offset by one so that zero
counts remain visible) and are annotated with the exact counts [TK: novel eGene,
sGene and pGene counts per biotype]. **(E)** Cumulative distribution of the
number of statistically independent signals per eGene in the combined
multi-ancestry expression fine-mapping. The *y* axis is the percentage of all
fine-mapped eGenes carrying at least *k* SuSiE (sum of single effects) credible
sets, out of every eGene with at least one credible set; the final point is the
percentage with five or more [TK: percentage of eGenes with ≥ 2 credible sets].
**(F)** Distribution of credible-set size — the number of variants contained in
a SuSiE credible set — for protein QTL (pQTL) fine-mapping in each ancestry
group and in the combined analysis (AoU), on a log10 scale. Groups are ordered
by increasing mean credible-set size. Boxes show the median and interquartile
range with whiskers extending to 1.5 times the interquartile range and
individual credible sets beyond that plotted as points. Black circles mark the
median credible-set size after restricting to credible sets with more than one
variant whose lead variants have an alternate-allele frequency of at least 1%.
The blue diamond marks the reported average credible-set size for UKBB pQTLs
[TK: 9.7 is hard-coded in `panels/panelF.R`; supply and cite the source of this
reference value] [TK: median credible-set size per group]. **(G)** Percentage of
credible sets resolved to a single putatively causal variant in the combined
multi-ancestry fine-mapping, by assay. A credible set was counted as resolved if
it contained a variant with a posterior inclusion probability (PIP) above 0.9;
the denominator is all credible sets for that assay. The number above each bar
is the count of resolved credible sets [TK: percentage and count of
single-variant credible sets for expression, protein and splicing]. **(H)**
Functional enrichment of putatively causal fine-mapped variants, by assay.
Within each assay's combined multi-ancestry SuSiE table, variants with
PIP > 0.9 were compared with a background of variants with PIP < 0.01 (each
variant represented once, at its maximum PIP across credible sets), and
enrichment for each annotation was tested with a two-sided Fisher exact test.
Annotations comprise Ensembl Variant Effect Predictor (VEP) consequences and
ENCODE candidate cis-regulatory element (cCRE) classes — promoter-like (PLS),
proximal enhancer-like (pELS) and distal enhancer-like (dELS) — grouped into
coding, regulatory and splice categories. Points give the odds ratio and bars
the 95% confidence interval from the Fisher exact test; the dashed line marks no
enrichment (odds ratio of 1) and the *y* axis is on a log10 scale. *P* values
were adjusted across all annotation-by-assay tests by the Benjamini–Hochberg
procedure and only annotations with an adjusted *P* below 0.05 are shown; the
CTCF-only cCRE class was excluded [TK: which annotations are displayed in the
rendered panel and their odds ratios].

Throughout, ancestry groups are AFR (African), AMR (admixed American), EAS (East
Asian), EUR (European), MID (Middle Eastern), SAS (South Asian) and, where
present, ADMIX (admixed/unassigned); AoU denotes the combined multi-ancestry
analysis of all groups together, and the same colour scheme is used in (B), (C)
and (F). Assays are expression (eQTL), splicing (sQTL) and protein abundance
(pQTL), and are coloured consistently in (D), (G) and (H). Error bars appear
only in (H), where they are 95% confidence intervals on the Fisher exact odds
ratio; no other panel plots an interval, and the boxplots in (F) summarise the
observed distribution rather than any uncertainty. Significance thresholds are
FDR < 0.05 in (C), Storey *q* < 0.05 in (D), PIP > 0.9 for putatively causal
variants in (G) and (H), and Benjamini–Hochberg adjusted *P* < 0.05 in (H). All
molecular phenotype and genotype data derive from the *All of Us* Research
Program controlled-access dataset [TK: data release / callset version] and were
analysed under its Data Use Agreement.

## Supplemental figure — cross-ancestry meta-analysis (MetaSoft)

**Figure S1. Cross-ancestry consistency of cis-QTL effects.** Per-ancestry
cis-QTL effect sizes and standard errors from the permutation results were
meta-analysed with MetaSoft in m-value mode, separately for expression, splicing
and protein abundance; the combined multi-ancestry analysis was excluded and
only phenotypes with a complete set of estimates across all six ancestry groups
(AFR, SAS, EUR, AMR, EAS, MID) were retained. **(A)** Proportion of genes, out
of all meta-analysed genes for that assay, by the number of ancestry groups
whose posterior probability that the effect is present (M value) is at least
0.9. Bars are grouped by assay. **(B)** Distribution of M values in each
ancestry group, faceted by assay; boxes show the median and interquartile range
with whiskers extending to 1.5 times the interquartile range [TK: number of
genes meta-analysed per assay; proportion of genes supported in all six groups].
This figure is produced by `Fig1/supp_metasoft.R` and is not part of the main
figure. [TK: confirm the supplemental figure number.]

## Open items

Placeholders that cannot be resolved from code (no access to the
controlled-access inputs; the panels are built from data absent from this
repository):

- `[TK: describe the externally supplied schematic]` — panel A is
  `blank_plot` in `panels/panelA.R`; the diagram is composed outside this repo.
- `[TK: n participants per ancestry group and the total]` (B) — counts come from
  the *All of Us* ancestry predictions intersected with the COMB eQTL sample
  list at run time.
- `[TK: number of phenotypes tested per assay; sample sizes for the protein and
  splicing assays]` (C) — see the sample-size caveat below.
- `[TK: confirm what the number printed beside the AoU point should represent]`
  (C) — `SummarizeSigQTLs()` in `00_setup.R` sets `num_sig = dplyr::n()`, i.e.
  the number of rows (phenotypes tested) in the group–assay stratum, **not** the
  number of significant QTLs the variable name implies. The AoU label in
  `panels/panelC.R` prints this value. Either the name or the computation is
  wrong; resolve before the legend claims either.
- `[TK: novel eGene, sGene and pGene counts per biotype]` (D) — `NOTES.md`
  records provisional counts from one staged run (protein-coding 1,207 / 2,035 /
  1,783; lncRNA 8,696 / 415 / 0; other non-coding RNA 246 / 0 / 0). These are
  not verifiable here and should be re-read from the rendered panel.
- `[TK: percentage of eGenes with ≥ 2 credible sets]` (E).
- `[TK: median credible-set size per group]` (F).
- `[TK: 9.7 is hard-coded ... supply and cite the source]` (F) — `ref_means` in
  `panels/panelF.R` fixes the UKBB pQTL comparison point at 9.7 with no
  derivation and no input file; `NOTES.md` confirms it is "still not derived
  from data". The legend must not present it as computed here.
- `[TK: percentage and count of single-variant credible sets ...]` (G).
- `[TK: which annotations are displayed in the rendered panel and their odds
  ratios]` (H) — the displayed set depends on which annotations survive the FDR
  filter in the data; the code whitelists twelve labels (splice donor, splice
  acceptor, splice region, stop gained, missense, inframe insertion, inframe
  deletion, 5' UTR, 3' UTR, promoter/PLS, proximal enhancer/pELS, distal
  enhancer/dELS) but any of them may be absent from the rendered figure.
- `[TK: data release / callset version]` — the manifest records bucket paths
  only; no release identifier appears in the code.
- `[TK: number of genes meta-analysed per assay; proportion ... in all six
  groups]` and `[TK: confirm the supplemental figure number]` (Fig. S1) — the
  supplemental script writes `xGeneMetasoftSummary.pdf` and carries no
  manuscript numbering.

Additional ambiguities and code/documentation conflicts:

- **Credible-set coverage level is never stated.** The SuSiE tables are read
  pre-computed; nothing in `Fig1/` sets or records the coverage (conventionally
  95%). Confirm from the fine-mapping pipeline before writing "95% credible set"
  in (E), (F) or (G).
- **Panel C uses one set of sample sizes for all three assays.** `panels/panelC.R`
  joins the per-group counts derived from the *expression* (COMB eQTL) sample
  list onto every assay, so the protein and splicing facets are plotted against
  expression-cohort sample sizes. The protein permutation filenames in
  `inputs.tsv` (`AoU_10750_*_pQTL.cis_qtl.txt.gz`) imply a different cohort size.
  Either this is intended (a shared genotyped cohort) or the *x* axis of the
  protein facet is mis-scaled; resolve before publication.
- **"Proportion of genes" in (C) is a proportion of phenotypes.** The
  denominator is rows of the permutation file, which for splicing are intron
  clusters rather than genes. The axis label in `PlotPropQTLSig()` says "genes".
- **Two different significance definitions within the figure.** Panel C
  recomputes Benjamini–Hochberg FDR from `pval_beta` in R; panel D uses the
  `qval` (Storey) column emitted by tensorQTL. Both are described as *q* < 0.05
  informally but they are not the same quantity.
- **`mean_cs` in panel F is a median.** `common_var_means` computes
  `median(cs_size)` into a column named `mean_cs`; the black points are medians.
  The boxplot ordering, separately, does use the mean.
- **README is stale on panel D.** `README.md` states the novel-pGene bar is
  "deferred pending a UKBB pQTL reference". The code builds it: `ukbb_pqtl` is in
  `inputs.tsv` and `panels/panelD.R` calls `LoadUKBBpQTL()`. `NOTES.md` records
  the resolution. The stale comment above `LoadUKBBpQTL()` in `00_setup.R`
  ("Currently deferred (no input file)") should also be removed.
- **README is stale on the panel H inputs.** `README.md` lists the `annot` key as
  `annotated/AoU_COMB_9k_{e,s}QTL_SusieMerged.annotated.tsv` plus
  `annotated/pQTL_susie_annotated.tsv`; `inputs.tsv` now points at
  `COMB.{expression,splicing,proteomics}_SusieMerged.annotated.ancestry_skew.tsv.gz`,
  shared with Figures 2 and 3.
- **README is stale on every staged path.** It describes `permutations/`,
  `finemapping/`, `reference/` and `external/` subdirectories; the manifest
  stages under `Fig1Data/` (e.g. `Fig1Data/eQTL/AFR.cis_qtl.txt.gz`), and
  `gene_gtf` / `ancestry_preds` resolve to bare filenames at the staging root,
  not under `reference/`.
- **`NOTES.md` is slightly stale on the panel H whitelist.** It describes
  `mature_mirna_variant` / `mature_miRNA_variant` as part of the twelve-label
  whitelist; the current `annotation_labels` in `PlotFinemappingEnrichment()`
  contains no microRNA entry, so no "Noncoding RNA" facet can be rendered even
  though `variant_functional_groups` defines that group.
