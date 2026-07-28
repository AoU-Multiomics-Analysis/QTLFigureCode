# Supplementary Figure 1 — draft legend

> Draft legend generated from the panel code. Values marked `[TK: ...]` must be
> filled in from the rendered figure or the data — they are not recoverable from
> the code alone. Verify against the final rendered figure before submission.

## Figure

**Supplementary Figure 1. Concordance of *All of Us* bulk RNA-seq gene expression
with GTEx whole blood across gene biotypes.**

Mean expression per gene was computed across all samples in each cohort as
transcripts per million (TPM), genes were matched between cohorts on Ensembl gene
identifier after stripping the version suffix, and gene biotype was assigned from
the GENCODE v45 annotation; only genes quantified in both cohorts are shown.
**(A)** Hexagonally binned density of protein-coding gene expression in *All of
Us* (AoU; x axis, log2(mean TPM + 1)) against GTEx v10 whole blood (y axis,
log2(mean TPM + 1)), restricted to protein-coding genes with a mean TPM greater
than 1 in at least one cohort. Hexagon fill encodes the log10 of the number of
genes falling in each bin; the red line is the identity line (y = x). Inset text
gives the Pearson correlation coefficient between the two log-transformed means
(*r* = [TK: correlation value printed in panel A]; n = [TK: number of
protein-coding genes plotted]). **(B)** Number of the protein-coding genes in
**A** that were detected in exactly one cohort, defined as a mean TPM greater
than 1 in one cohort and below 1 in the other; bar height is the count of unique
genes specific to *All of Us* and to GTEx, respectively ([TK: cohort-specific
gene counts]). **(C)** Gene Ontology (GO) enrichment of the two cohort-specific
protein-coding gene sets defined in **B**, computed with Enrichr against the GO
Molecular Function, GO Cellular Component and GO Biological Process 2025 term
libraries, which were pooled before ranking. The nine most significant terms per
cohort, ordered by Enrichr's adjusted *P* value, are shown; facets separate
*All of Us*-specific from GTEx-specific gene sets. The x axis is the negative
log10 of the adjusted *P* value (false-discovery rate, FDR) and the vertical
reference line is drawn at 1.3, corresponding to an FDR of approximately 0.05.
**(D)** Expression of transcription-factor (TF) genes in the two cohorts. Genes
were retained if they were protein-coding, carried a symbol present in the Human
Transcription Factors list (v1.01), and had a mean TPM greater than 1 in at least
one cohort ([TK: number of TF genes retained]); the y axis is log2(mean TPM + 1)
and each box shows the median and interquartile range with whiskers extending to
the most extreme value within 1.5× the interquartile range (points beyond are
plotted individually). Expression was compared between cohorts with a two-sided
paired *t*-test on the log2(mean TPM + 1) values (*P* = [TK: paired t-test P
value]) and by ordinary least-squares regression of GTEx on *All of Us*
log2(mean TPM + 1) (slope = [TK], *P* = [TK]); both statistics are written to the
run log and are not annotated on the panel. **(E)** As in **A**, for long
non-coding RNA (lncRNA) genes, restricted to lncRNAs with a mean TPM greater than
0.1 in at least one cohort (*r* = [TK: correlation value printed in panel E]; n =
[TK: number of lncRNA genes plotted]). **(F)** As in **B**, for the lncRNA genes
in **E**: the number of lncRNAs with a mean TPM greater than 1 in one cohort and
below 1 in the other ([TK: cohort-specific lncRNA counts]).

Throughout, bar heights in **B** and **F** are gene counts rather than estimated
quantities, and no error bars or confidence intervals are drawn on any panel.
Cohort is colour-coded consistently in **B**, **D** and **F** using the
CARTOColors "Safe" qualitative palette; the hexagon fill in **A** and **E** uses
the viridis continuous scale. Expression data are *All of Us* bulk RNA-seq TPMs
(controlled-access; [TK: number of *All of Us* participants/samples and the
sequenced tissue]) and GTEx v10 whole-blood TPMs ([TK: number of GTEx samples]),
annotated with GENCODE v45. Abbreviations: TPM, transcripts per million; AoU,
*All of Us*; GTEx, Genotype-Tissue Expression project; lncRNA, long non-coding
RNA; TF, transcription factor; GO, Gene Ontology; FDR, false-discovery rate;
IQR, interquartile range.

## Open items

- `[TK: correlation value printed in panel A]` and `[TK: ... panel E]` — the
  Pearson coefficient is computed at run time by `cor()` in
  `plot_gene_comparison()` (`SuppFig1/00_setup.R:103-105`) and drawn as an
  `annotate("text")` label; the number exists only in the rendered figure.
- `[TK: number of protein-coding genes plotted]`, `[TK: number of lncRNA genes
  plotted]` — set by the panel filters (`SuppFig1/panels/panelA.R:5-7`,
  `SuppFig1/panels/panelD.R:5-7`) applied to controlled-access data.
- `[TK: cohort-specific gene counts]`, `[TK: cohort-specific lncRNA counts]` —
  the bar heights are `dplyr::count()` outputs computed at run time
  (`SuppFig1/00_setup.R:125-133`).
- `[TK: number of TF genes retained]` — the row count of
  `tf_expression_table()` (`SuppFig1/00_setup.R:212-217`).
- `[TK: paired t-test P value]`, `[TK: regression slope]`, `[TK: regression P]` —
  `SuppFig1/panels/panelC.R:8-13` computes a `broom::tidy(lm(...))` and a paired
  `t.test()` and only `message()`/`print()`s them; they are not rendered on the
  panel, so if the manuscript quotes them they must be read off a pipeline run.
  Note also that `t.test()` is called with default arguments, i.e. two-sided.
- `[TK: number of *All of Us* participants/samples and the sequenced tissue]` —
  the manifest gives only the file
  (`expression/aou_rnaseq_121624.gene_tpm.fixed.gct.gz`,
  `SuppFig1/inputs.tsv:2`); the sample count is the column count of a
  controlled-access matrix. The AoU tissue is never stated in code — only
  `SuppFig1/README.md` implies whole blood by comparing to GTEx whole blood.
- `[TK: number of GTEx samples]` — column count of the public GTEx v10
  whole-blood GCT; not recorded anywhere in the repo.
- Enrichr's statistical test and multiple-testing correction are **not** set in
  this repo — `SuppFig1/00_setup.R:189-195` only calls `enrichR::enrichr()` and
  the panel consumes the returned `Adjusted.P.value`. If the manuscript is to
  name the test (Fisher's exact) or the correction (Benjamini–Hochberg), cite
  Enrichr rather than this code. The `-log10(FDR)` axis label
  (`SuppFig1/panels/panelB.R:13`) is the code's own wording for that column.
- Which of the three GO ontologies each displayed term came from is not
  recoverable from the panel: the three libraries are `bind_rows()`-ed before the
  top-9 slice (`SuppFig1/00_setup.R:192-194`,
  `SuppFig1/panels/panelB.R:8`), and the ontology is not shown on the axis.
- **Panel-letter collision between the docs and the rendered figure.** The panel
  *files* are A–D, but `assemble.R` letters six leaf plots, so the rendered tags
  are: file A → **A, B**; file B (GO) → **C**; file C (TF) → **D**; file D
  (lncRNA) → **E, F** (`SuppFig1/assemble.R:12-21`, pinned by
  `tests/test_suppfig1_assembly.R`). `SuppFig1/README.md:49` and every heading in
  `SuppFig1/NOTES.md` call the TF panel "Panel C" and the GO work "Panel B",
  which are the *file* letters and contradict the printed tags. This legend uses
  the rendered tags. The README's own summary table (`README.md:15-20`) uses the
  rendered tags and is correct.
- **Detection threshold is not matched to the inclusion threshold in E/F.**
  `panelD.R` pre-filters lncRNAs at mean TPM > 0.1 but calls
  `plot_gene_comparison()` without `detection_threshold`, so the bars in **F**
  still use the default threshold of 1 (`SuppFig1/00_setup.R:100`,
  `SuppFig1/panels/panelD.R:5-8`). Panels A/B use 1 for both. Confirm this is
  intended before describing **F**; the design spec
  (`docs/superpowers/specs/2026-07-24-suppfig1-cleanup-design.md:128`) says only
  "detection counts" and does not resolve it.
- **Panel D's x axis is unlabelled in code.** `panelC.R:19` maps `x = name` and
  sets only `ylab()`, so the rendered x-axis title is the literal column name
  `name` rather than "Cohort" (the detection panels do set `xlab("Cohort")`,
  `SuppFig1/00_setup.R:132`). Cosmetic, but it will appear in the figure.
- **Naming.** The repo, the output file (`SupplementalFig1.pdf`,
  `SuppFig1/assemble.R:23`) and all docs say "Supplemental Figure 1"; this legend
  follows the manuscript's "Supplementary Figure 1". Align before submission.
- **Spec-vs-code.** The design spec lists `merged_tpms()` as feeding "all extras"
  (`...-suppfig1-cleanup-design.md:107`), but
  `SuppFig1/extra/gene_type_detection.R` deliberately reads the AoU GCT directly
  so as not to narrow the gene universe to the AoU∩GTEx intersection
  (documented in `SuppFig1/NOTES.md:31-37`). This does not affect the assembled
  figure — none of the `extra/` scripts (snRNA/miRNA comparisons, rRNA histogram,
  hemoglobin-gene scatter, per-biotype detection counts) appear in
  `SupplementalFig1.pdf`, and none are described in this legend.
- `go_annotations` is declared in `SuppFig1/inputs.tsv:6` but consumed by nothing
  in the pipeline, so no panel derives from it despite it being controlled-access.
- This figure contains **no** ancestry-group or QTL stratification, so none of the
  manuscript's population labels (AFR, AMR, EAS, EUR, MID, SAS) or QTL
  abbreviations (eQTL, pQTL, sQTL) apply here. If the manuscript's supplementary
  legends share a boilerplate abbreviation block, trim it for this figure.
