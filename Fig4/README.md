# Figure 4 — QTL–GWAS colocalization and rare-variant TWAS

Code to generate **Figure 4** for the manuscript. The figure connects molecular
QTL colocalization with GWAS traits to downstream rare-variant signal: how many
genes colocalize per trait category and at what lenient credible-set rate, the R²
that rare variants add to TWAS, whether colocalized genes carry inflated
conditional TWAS signal versus matched controls, and whether that signal enriches
in clinically important gene sets.

> **Data note.** The inputs derive from *All of Us* **controlled-access** data
> and are **not** in this repository. Run this on the Researcher Workbench.

## The figure

| Panel | Content | Source |
|-------|---------|--------|
| **A** | Colocalizing-gene counts + **lenient** credible-set coloc rate, by trait category (traits with > 5 colocalizing genes; top 10 per category; two categories dropped) | `harmonized_coloc.gene.tsv.gz` + `coloc_rate_by_trait.tsv` |
| **B** | Histogram of the R² increase (%) contributed by rare variants (x truncated to −20…40 %) | `R2_TWAS_combined.tsv` |
| **C** | Common/rare conditional-TWAS QQ, colocalized vs. matched non-colocalized (blood traits) | `ConditionalTWAS.cleaned.tsv.gz` |
| **D** | Odds-ratio forest: rare/common TWAS-signal enrichment across gene sets | conditional TWAS + ClinGen + drug interactions |

### Panel A — which rate, and what is filtered out

The lollipop is the **lenient** union credible-set coloc rate
(`pct_lenient_union`): `coloc_by_trait()` in `00_setup.R` calls
`build_plot_data(..., rate_type = "lenient")`. Note that `build_plot_data()`'s own
default is `rate_type = "stringent"` (`pct_stringent_union`) — anything else that
calls the helper directly gets the *stringent* rate unless it asks for lenient.

`panels/panelA.R` then trims the table before plotting: it drops the
`neuropsychiatric` and `renal_urological` categories, keeps **traits** with
`count > 5` colocalizing genes, and keeps the top ten **traits within each
category** by count. The "top 10" is therefore per category, not a global top-10
of categories.

### Panel B — truncated axis

`panels/panelB.R` sets `xlim(-20, 40)` on the histogram. That is a scale limit,
not a zoom: genes whose R² increase falls outside −20…40 % are **removed** from
the data before binning rather than piling up at the edge, so the plotted gene
count is smaller than the number of genes passing `Pval < 1e-5`.

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

### Notebook

`Fig4/Fig4.ipynb` is an annotated R-kernel notebook that reproduces this figure
step by step, with markdown explaining each analysis. It is generated from the
pipeline scripts (so it stays faithful) and regenerated with:

    python3 tools/build_notebooks.py Fig4

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
