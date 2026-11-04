# QTLFigureCode

The [data preparation guide](preparation/README.md) documents the filtering,
identifier matching, counts, statistical summaries, and scripts that produce
the plot inputs. 

R code that generates the figures for a multi-omic **molecular-QTL** manuscript:
cis-QTL discovery, fine-mapping, and cross-assay colocalization of **expression
(eQTL)**, **protein (pQTL)** and **splicing (sQTL)** QTLs across genetic-ancestry
groups in the [*All of Us*](https://allofus.nih.gov/) cohort.


> [!IMPORTANT]
> **The inputs are *All of Us* controlled-access genomic data and are not in
> this repository.**  Access is via the *All of Us* Researcher Workbench under its Data Use Agreement. See
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
└── README.md               # (this file)
```

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
staging root that defaults to the working directory), plus
`source_repo` / `source_workflow` recording which upstream
repository and workflow produced the file (`-` where unattributed, `external`
for third-party references). The public manifests use `-` in `gs_uri` for
private GCS locations. Public HTTPS sources remain available. Supply private
inputs at their `local_subpath` before running the code. See the
[private manifest guide](docs/private-input-manifests.md) to use a separate local
manifest with the original addresses.

Any other column is treated as metadata and attached to the rows a load returns.
When a separate manifest supplies an address, a file is fetched **only when the local copy is missing** — via
`gsutil` for `gs://`, or `download.file` for `http(s)://` and `file://` — so
pre-staged inputs skip the copy. Panels reach inputs through `input_path(key)`,
`load_input(key)` and `input_files(key)`, and each panel declares only what it
needs. To add or re-point an input you edit `inputs.tsv` — no script changes.

**Notebooks.** `tools/build_notebooks.py` generates the current notebooks from
`notebooks/current/` and expands their R helper includes into visible cells.
They show source-result processing and plot construction. Regenerate all with
`python3 tools/build_notebooks.py`, or one with
`python3 tools/build_notebooks.py Fig4`. Run them again after a code change.
Use `--legacy` for the older modular-pipeline notebooks. The
[current notebook guide](notebooks/current/README.md) gives input requirements
and execution commands.

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



```bash
Rscript tests/run_all.R
```

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
- **Python 3.9 or later** for `tools/build_notebooks.py` (standard library only).
- External tools where noted: **MetaSoft** (Java) for Fig1's supplemental;
  **FastENLOC** produces Fig3's inputs upstream (Fig3 treats them as inputs, not
  as something it runs).

---

## Data & inputs

Access to the underlying data is governed by the *All of Us* Data Use Agreement
and is **not** provided through this repository. The pipelines here have been
exercised against unit-test mock fixtures; the real end-to-end run happens on the
Workbench, where the data resides.


