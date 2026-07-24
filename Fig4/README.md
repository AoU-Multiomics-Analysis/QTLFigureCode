# Figure 4 — QTL–GWAS colocalization and rare-variant TWAS

Code to generate **Figure 4** for the manuscript. The figure connects molecular
QTL colocalization with GWAS traits to downstream rare-variant signal: how many
genes colocalize per trait category and at what credible-set rate, the R² that
rare variants add to TWAS, whether colocalized genes carry inflated conditional
TWAS signal versus matched controls, and whether that signal enriches in
clinically important gene sets.

> **Data note.** The inputs derive from *All of Us* **controlled-access** data
> and are **not** in this repository. Run this on the Researcher Workbench.

## The figure

| Panel | Content | Source |
|-------|---------|--------|
| **A** | Colocalizing-gene counts + credible-set coloc rate, by trait category | `harmonized_coloc.gene.tsv.gz` + `coloc_rate_by_trait.tsv` |
| **B** | Histogram of the R² increase (%) contributed by rare variants | `R2_TWAS_combined.tsv` |
| **C** | Common/rare conditional-TWAS QQ, colocalized vs. matched non-colocalized (blood traits) | `ConditionalTWAS.cleaned.tsv.gz` |
| **D** | Odds-ratio forest: rare/common TWAS-signal enrichment across gene sets | conditional TWAS + ClinGen + drug interactions |

## Repository layout, running, tests

Run from the **project root**:

    Rscript fig4/run_all.R          # panels -> assemble -> Figure4.pdf / .png
    Rscript fig4/panels/panelD.R    # rebuild one panel
    Rscript fig4/assemble.R         # recompose
    rm -rf fig4/data && Rscript fig4/run_all.R   # force-rebuild the derived cache
    Rscript tests/run_all.R         # unit tests (no controlled-access data)

The panels share upstream work through memoized builders in `00_setup.R`
(`read_conditional_twas`, `gene_sets`, `coloc_by_trait`, `qq_match_source`,
`gene_twas_data`), cached to `fig4/data/` **by name** — clear that directory
after changing an input.

## Requirements

R ≥ 4.3, and CRAN packages `tidyverse`, `data.table`, `patchwork`, `ggh4x`,
`MatchIt`, `arrow`, `broom`, `scales`, `showtext`, `sysfonts`, `paletteer`.
(`readxl` is **not** needed — the GeneBayes read was unused and removed.)

## Data & inputs

Every input is declared in `fig4/inputs.tsv`. The two `high_level_coloc_outputs/`
tables come from a single archive:

    gsutil cp gs://fc-secure-…/high_level_coloc_outputs.tar.gz .
    tar -xf high_level_coloc_outputs.tar.gz

The ClinGen `reference/evidence_clingen/` parquet, `reference/interactions.tsv`,
`R2_TWAS_combined.tsv` and `ConditionalTWAS.cleaned.tsv.gz` are local-only and
must be pre-staged. See [`NOTES.md`](NOTES.md) — in particular the
result-affecting fix to the duplicated conditional-TWAS table (Panel D).

## Provenance

This directory replaces a single 774-line notebook export (`Fig4.R`, recoverable
from git history). The restructuring and its fixes are recorded in
[`docs/superpowers/specs/2026-07-24-fig4-cleanup-design.md`](../docs/superpowers/specs/2026-07-24-fig4-cleanup-design.md).
