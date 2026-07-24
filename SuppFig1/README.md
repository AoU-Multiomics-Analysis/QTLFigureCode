# Supplemental Figure 1 — All of Us vs GTEx expression

Code to generate **Supplemental Figure 1**: a comparison of *All of Us* bulk
RNA-seq expression to GTEx whole blood across gene types.

> **Data note.** `aou_tpms` (and `go_annotations`) derive from *All of Us*
> **controlled-access** data and are not in this repository. The GTEx TPMs,
> GENCODE GTF and TF list are public and self-stage from `https://` URIs.

## The figure

`SupplementalFig1.pdf` — a 2×2 of composites, lettered A–F (panels A and D are
each a 2-plot pair):

| Tags | Content |
|------|---------|
| **A, B** | protein-coding AoU-vs-GTEx TPM correlation (hex) + cohort-specific detection counts |
| **C** | GO enrichment of cohort-specific genes (AoU vs GTEx) |
| **D** | transcription-factor expression (AoU vs GTEx), log2(TPM+1) |
| **E, F** | lncRNA AoU-vs-GTEx TPM correlation (hex) + detection counts |

Additional per-gene-type comparisons (snRNA, miRNA, rRNA, hemoglobin genes, and
per-type detection counts) are in `extra/`, each saving its own PDF.

## Running

From the **project root**:

    Rscript suppfig1/run_all.R              # panels -> assemble -> SupplementalFig1.pdf/.png
    Rscript suppfig1/extra/rrna_expression.R  # a standalone extra
    rm -rf suppfig1/data && Rscript suppfig1/run_all.R   # force-rebuild the cache
    Rscript tests/run_all.R                 # unit tests (no controlled-access data)

The panels share upstream through memoized builders in `00_setup.R`
(`gene_annotations`, `merged_tpms`, `tf_list`, `go_enrichment`), cached to
`suppfig1/data/` **by name** — clear it after changing an input. `go_enrichment`
makes a live enrichR HTTP call, cached so it runs once.

## Requirements

R ≥ 4.3, and CRAN/Bioconductor packages `tidyverse`, `data.table`, `rtracklayer`,
`enrichR`, `hexbin`, `ggrepel`, `rcartocolor`, `broom`, `patchwork`, `showtext`,
`sysfonts`, `paletteer`.

## Data & inputs

Every input is declared in `suppfig1/inputs.tsv`; the three public inputs
self-stage over https. See [`NOTES.md`](NOTES.md), in particular the
result-affecting TF double-log fix (Panel C).

## Provenance

Replaces a single 316-line notebook export (`SuppFig1.R`, recoverable from git
history). The restructuring and its fixes are recorded in
[`docs/superpowers/specs/2026-07-24-suppfig1-cleanup-design.md`](../docs/superpowers/specs/2026-07-24-suppfig1-cleanup-design.md).
