# QTLFigureCode

The [current notebook guide](notebooks/current/README.md) describes the November 2,
2026 main figures and the 16 supplement notebooks. It lists local input paths,
execution steps, source-data coverage, and the limits of each reconstruction.
The [data preparation guide](preparation/README.md) documents the filtering,
identifier matching, counts, statistical summaries, and scripts that produce
the plot inputs. Source data remain outside Git.

R code that generates the figures for a multi-omic **molecular-QTL** manuscript:
cis-QTL discovery, fine-mapping, and cross-assay colocalization of **expression
(eQTL)**, **protein (pQTL)** and **splicing (sQTL)** QTLs across genetic-ancestry
groups in the [*All of Us*](https://allofus.nih.gov/) cohort.

Each manuscript figure has its own directory with a self-contained pipeline. The
analyses run on the *All of Us* Researcher Workbench (Terra/Jupyter), where the
input data lives.

> [!IMPORTANT]
> **The inputs are *All of Us* controlled-access genomic data and are not in
> this repository.** Nothing derived from a participant — staged input files,
> cached intermediates, or rendered figures — may be committed. Access is via
> the *All of Us* Researcher Workbench under its Data Use Agreement. See
> [Data & inputs](#data--inputs).

---

## The figures

All six figures follow the same modular, manifest-driven template.

| Dir | Figure | What it shows | Panels |
|-----|--------|---------------|--------|
| [`Fig1/`](Fig1/README.md) | **1** | cis-QTL discovery across the three assays and ancestry groups: cohort composition, proportion of genes with a significant QTL by sample size, novelty vs. prior consortia, and SuSiE fine-mapping architecture | A–H |
| [`Fig2/`](Fig2/README.md) | **2** | properties of fine-mapped causal variants (SuSiE, PIP > 0.9): allele-frequency spectrum, ancestry enrichment, molecular effect sizes, and links to expression-outlier / dosage phenotypes | A–G |
| [`Fig3/`](Fig3/README.md) | **3** | how often a signal is **shared between molecular layers** — FastENLOC colocalization between each assay pair at gene, credible-set and variant level, and whether sharing depends on allele frequency | A–E |
| [`Fig4/`](Fig4/README.md) | **4** | colocalization of molecular QTLs with GWAS traits, grouped by trait category, plus rare-variant TWAS | A–D |
| [`Fig5/`](Fig5/README.md) | **5** | rare-variant eQTL constraint and sequence-model variant-effect dissection: fine-mapped rare ClinGen variants, plus worked *ACVRL1* (PromoterAI) and *SPR* (Borzoi) examples | A–G |
| [`SuppFig1/`](SuppFig1/README.md) | **S1** | *All of Us* vs GTEx whole-blood expression: cohort-specific genes, GO enrichment, TF expression | A–D |

Each figure directory has its own `README.md` and `NOTES.md` — the latter
records deferred analyses and every way the cleaned code deliberately differs
from the original. **Read a figure's `NOTES.md` before comparing its output to a
previously published version**; several changes affect results (notably Fig3).

Some figures also ship standalone supplemental scripts (`supp_*.R`) that are not
part of `run_all.R`: Fig1's MetaSoft figure, Fig2's ancestry-composition figure,
and Fig3's colocalization-QC and per-gene-proportion figures. Fig5 keeps a
larger set of retained provenance and side analyses under
[`Fig5/extra/`](Fig5/extra/README.md).

---

## Repository layout

```
.
├── common/
│   └── manifest.R          # shared input-manifest + staging (used by every figure)
├── Fig1/ … Fig5/           # one directory per manuscript figure
├── SuppFig1/               # supplemental figure, same template
├── tests/                  # unit + regression tests, runnable without any real data
│   ├── fixtures.R          #   mock project builders (Fig3, Fig4, Fig5, SuppFig1)
│   ├── helpers.R           #   ok() / expect_error()
│   ├── run_all.R           #   runs every tests/test_*.R
│   └── test_*.R
├── tools/
│   └── build_notebooks.py  # regenerates every Fig<N>.ipynb from the pipeline scripts
├── docs/
│   ├── figure-cleanup-playbook.md   # the workflow for cleaning up a figure
│   └── superpowers/                 # design specs + implementation plans per cleanup
├── CLAUDE.md               # guidance for AI sessions working in this repo
└── README.md               # (this file)
```

Staged inputs, per-panel object caches, and rendered figures all land in
git-ignored locations (see [`.gitignore`](.gitignore) and each figure's own
`.gitignore`); they are never committed. The one exception is
`Fig5/reference/` — a handful of small **public** reference files (ClinGen gene
validity, ENCODE cCRE regions, gene models) that are committed deliberately so
Fig5's worked examples build without staging.

---

## Shared conventions

Every figure follows one template — Fig2 established the shape, and the manifest
system was layered on across the later cleanups:

- **`00_setup.R`** — libraries, env-overridable paths, the `showtext` font and a
  single `DPI` constant, a `theme_set`/`theme_update` design system with `SZ_*`
  size constants, and any analysis/plot helpers folded in so the figure is
  self-contained. Figures with a heavy shared upstream (Fig3, Fig4, Fig5) also
  define memoized builders here that cache to `data/`.
- **`inputs.tsv`** — the input manifest (see below).
- **`panels/panel<A…>.R`** — each panel sources setup, loads only the inputs it
  needs, builds one `ggplot`/`patchwork` object, and caches it to `objects/`.
- **`assemble.R`** — composes the cached panel objects (via
  [patchwork](https://patchwork.data-imaginist.com/)) and writes the final
  `Figure<N>.pdf` plus a PNG preview. **Panel tagging happens here and only
  here.**
- **`run_all.R`** — panels → assemble, end to end.
- **`Fig<N>.ipynb`** — a generated notebook (below).
- **`NOTES.md`** — behaviour changes and deferred work.

**Manifest-driven inputs** (`common/manifest.R`). Each figure's `inputs.tsv`
maps a logical `key` to one or more files: a `local_subpath` (relative to a
staging root that defaults to the working directory) and, where known, a `gs://`
origin, plus `source_repo` / `source_workflow` recording which upstream
repository and workflow produced the file (`-` where unattributed, `external`
for third-party references). Any other column is treated as metadata and
attached to the rows a load returns. A file is fetched **only when the local copy is missing** — via
`gsutil` for `gs://`, or `download.file` for `http(s)://` and `file://` — so
pre-staged inputs skip the copy. Panels reach inputs through `input_path(key)`,
`load_input(key)` and `input_files(key)`, and each panel declares only what it
needs. To add or re-point an input you edit `inputs.tsv` — no script changes.

**Notebooks.** `tools/build_notebooks.py` generates one Jupyter notebook per
figure by pulling the code cells straight from the pipeline scripts, so the
notebooks cannot drift from the code. Regenerate them all with `python3
tools/build_notebooks.py`, or one with `python3 tools/build_notebooks.py Fig4`.
The tracked notebooks must round-trip through the generator — regenerating
should produce no diff. They are meant to run on the Workbench, where the inputs
live; the generator never executes them.

**Rendering.** Text uses the **Inter** font, fetched from Google Fonts at runtime
via [showtext](https://cran.r-project.org/package=showtext), with an automatic
fallback to the default sans font when Google Fonts is unreachable. Every output
DPI is matched to the `showtext` DPI so text scales correctly. (The value itself
is per figure: 300 for Fig1, Fig2 and SuppFig1; 400 for Fig3–Fig5.)

> [!NOTE]
> **Run every script from the project root**, e.g. `Rscript fig3/run_all.R`.
> Scripts refer to themselves with a lowercase `fig<N>/…` path even though the
> directories are capitalised (`Fig3/`). This works on the case-insensitive
> filesystems the code was developed on (macOS) and on the Workbench, but **would
> fail on a case-sensitive host** (most Linux) unless the path case is matched.

---

## Running

Pick a figure directory and run its pipeline from the project root. For example:

```bash
# Build every panel (staging inputs as needed) and assemble the figure
Rscript fig3/run_all.R

# Rebuild a single panel after editing it, then recompose
Rscript fig3/panels/panelD.R
Rscript fig3/assemble.R

# Supplemental figures are standalone
Rscript fig3/supp_coloc_qc.R
```

Fig2 is the one variant: it reads its raw inputs once in `fig2/01_prep_data.R`
and writes a prepared frame per panel, so the prep step must run before the
panels (`fig2/run_all.R` does both).

See each figure's `README.md` for its exact panel set, inputs, outputs, and any
figure-specific steps (Fig1's MetaSoft supplemental, Fig3–Fig5's derived-data
caches, etc.).

> [!TIP]
> The memoized builders in Fig3–Fig5 do **not** cascade `refresh = TRUE` through
> their dependencies. After changing an input or a builder, clear the whole
> cache (`rm -rf fig<N>/data`) rather than refreshing a single builder.

### Tests

`tests/` holds unit and regression checks that need **no controlled-access
data** — they run against mock projects built by `tests/fixtures.R`:

```bash
Rscript tests/run_all.R
```

The suite covers the shared manifest loader and the Fig3, Fig4, Fig5 and
SuppFig1 pipelines — including regression tests for defects found in the
original code and, for each of those figures, a test asserting the **rendered**
panel-letter sequence of the assembled figure. Fig1 and Fig2 have no coverage
yet; their prep steps are not factored into separately testable seams.

---

## Requirements

- **R** ≥ 4.3
- Core CRAN packages used across figures: `tidyverse`, `data.table`, `patchwork`,
  `scales`, `showtext`, `sysfonts`, `paletteer`.
- Additional per figure — e.g. `igraph` (Fig3); `broom` (Fig3, Fig4, SuppFig1);
  `readxl` (Fig1, Fig5); `ggrepel` (Fig1, Fig5, SuppFig1); `ggbeeswarm`, `ggdist`
  and `ggh4x` (Fig2); `cowplot` and `ggseqlogo` (Fig5); and the Bioconductor
  package `rtracklayer` (Fig1, Fig5, SuppFig1). Each figure's `README.md` lists
  its own dependencies and an install snippet.
- **`gsutil`** (Google Cloud SDK) for on-demand input staging from the secure
  bucket — present by default on the *All of Us* Workbench.
- **Python 3** for `tools/build_notebooks.py` (standard library only).
- External tools where noted: **MetaSoft** (Java) for Fig1's supplemental;
  **FastENLOC** produces Fig3's inputs upstream (Fig3 treats them as inputs, not
  as something it runs).

---

## Data & inputs

Every input derives from *All of Us* controlled-access data staged from a secure
Google Cloud Storage bucket on the Researcher Workbench. Each figure's
`inputs.tsv` is the record of its inputs' origins — but **not every input has a
`gs://` URI**. Those that don't must be pre-staged by hand; the pipeline cannot
self-stage them. Fig2, Fig4 and Fig5 currently have no `gs://` origins at all,
and Fig1, Fig3 and SuppFig1 have a mix. Each figure's `README.md` and `NOTES.md`
flag which of its inputs still need one.

Access to the underlying data is governed by the *All of Us* Data Use Agreement
and is **not** provided through this repository. The pipelines here have been
exercised against unit-test mock fixtures; the real end-to-end run happens on the
Workbench, where the data resides.

---

## Development notes

Every figure was refactored from a single flat notebook export into the modular
pipeline above. The workflow for doing another one is written up in
[`docs/figure-cleanup-playbook.md`](docs/figure-cleanup-playbook.md), with the
per-cleanup design specs and implementation plans under
[`docs/superpowers/`](docs/superpowers/). [`CLAUDE.md`](CLAUDE.md) points AI
sessions at the same material.

Each cleaned figure's `NOTES.md` records deferred work and any deliberate
difference from the original analysis — including several changes that affect
figure output, so read the relevant one before comparing against a previously
published version.
