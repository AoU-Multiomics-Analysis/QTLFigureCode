# Supplemental Figure 1 — All of Us vs GTEx expression

Code to generate **Supplemental Figure 1**: a comparison of *All of Us* bulk
RNA-seq expression to GTEx whole blood across gene types.

> **Data note.** `aou_tpms` (and `go_annotations`) derive from *All of Us*
> **controlled-access** data and are not in this repository. The GTEx TPMs,
> GENCODE GTF and TF list are public and self-stage from `https://` URIs.

## The figure

`SupplementalFig1.pdf` — a 2×2 of composites, lettered A–F (the `panelA.R` and
`panelD.R` files each contribute a 2-plot pair):

| Tags | Content |
|------|---------|
| **A, B** | protein-coding AoU-vs-GTEx TPM correlation (hex) + cohort-specific detection counts |
| **C** | GO enrichment of cohort-specific genes (AoU vs GTEx) |
| **D** | transcription-factor expression (AoU vs GTEx), log2(TPM+1) |
| **E, F** | lncRNA AoU-vs-GTEx TPM correlation (hex) + detection counts |

> **Panel *files* vs rendered *tags*.** `assemble.R` composes four panel objects,
> built by `panels/panelA.R`–`panels/panelD.R`, but
> `plot_annotation(tag_levels = "A")` letters *leaf* plots — and `panelA.R` and
> `panelD.R` each build a two-plot patchwork. So four files produce six tags:
> `panelA.R` → **A, B**; `panelB.R` (GO enrichment) → **C**; `panelC.R`
> (transcription factors) → **D**; `panelD.R` (lncRNA) → **E, F**. The six-tag
> layout is intentional and is pinned by `tests/test_suppfig1_assembly.R`.
> Throughout this README and [`NOTES.md`](NOTES.md), "Panel *X*" means the
> rendered tag, not the file.

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
makes a live enrichR HTTP call, cached so it runs once. That call only submits
the cohort-specific gene lists: the enrichment test and its multiple-testing
correction are Enrichr's, configured nowhere in this repo, and the returned
statistics (including `Adjusted.P.value`, which Panel C plots as `-log10(FDR)`)
are consumed as-is. Enrichr is a live service rather than a declared input, so it
has no row in `inputs.tsv`.

## Requirements

R ≥ 4.3, and CRAN/Bioconductor packages `tidyverse`, `data.table`, `rtracklayer`,
`enrichR`, `hexbin`, `ggrepel`, `rcartocolor`, `broom`, `patchwork`, `showtext`,
`sysfonts`, `paletteer`.

## Data & inputs

Every input is declared in `suppfig1/inputs.tsv`; the three public inputs
self-stage over https. Private GCS addresses are omitted. Supply `aou_tpms`
and `go_annotations` at their listed local paths, or use a separate local
manifest as described in the [private manifest guide](../docs/private-input-manifests.md).
See [`NOTES.md`](NOTES.md), in particular the
result-affecting TF double-log fix (rendered Panel D, built by
`panels/panelC.R`).

## Provenance

Replaces a single 316-line notebook export (`SuppFig1.R`, recoverable from git
history). The restructuring and its fixes are recorded in
[`docs/superpowers/specs/2026-07-24-suppfig1-cleanup-design.md`](../docs/superpowers/specs/2026-07-24-suppfig1-cleanup-design.md).
