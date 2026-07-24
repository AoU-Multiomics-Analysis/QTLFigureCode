# QTLFigureCode

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

| Dir | Figure | What it shows | Status |
|-----|--------|---------------|--------|
| [`Fig1/`](Fig1/README.md) | **1** | cis-QTL discovery across the three assays and ancestry groups: cohort composition, proportion of genes with a significant QTL by sample size, novelty vs. prior consortia, and SuSiE fine-mapping architecture | Modular pipeline + manifest |
| [`Fig2/`](Fig2/README.md) | **2** | properties of fine-mapped causal variants (SuSiE, PIP > 0.9): allele-frequency spectrum, ancestry enrichment, molecular effect sizes, and links to expression-outlier / dosage phenotypes | Modular pipeline (pre-manifest) |
| [`Fig3/`](Fig3/README.md) | **3** | how often a signal is **shared between molecular layers** — FastENLOC colocalization between each assay pair at gene, credible-set and variant level, and whether sharing depends on allele frequency | Modular pipeline + manifest |
| `Fig4/` | **4** | colocalization of molecular QTLs with GWAS traits, grouped by trait category | Flat script — not yet cleaned |
| `Fig5/` | **5** | — | Empty (planned) |

Each figure directory has its own `README.md` (and, where relevant, a `NOTES.md`
recording deferred analyses and any way the cleaned code deliberately differs
from the original). Start there for the detail on a given figure.

---

## Repository layout

```
.
├── common/
│   └── manifest.R          # shared input-manifest + gsutil staging (Fig1, Fig3, …)
├── Fig1/ … Fig5/           # one directory per manuscript figure
├── tests/                  # unit + regression tests, runnable without any real data
│   ├── fixtures.R          #   mock project builder used by the Fig3 tests
│   ├── helpers.R           #   ok() / expect_error()
│   ├── run_all.R           #   runs every tests/test_*.R
│   └── test_*.R
├── docs/superpowers/       # design specs and implementation plans for the cleanups
└── README.md               # (this file)
```

Staged inputs, per-panel object caches, and rendered figures all land in
git-ignored locations (see [`.gitignore`](.gitignore) and each figure's own
`.gitignore`); they are never committed.

---

## Shared conventions

Figures that have been cleaned up (**Fig1** and **Fig3**) follow one template —
Fig2 established its shape, and Fig1/Fig3 extended it with the manifest system:

- **`00_setup.R`** — libraries, env-overridable paths, the `showtext` font and a
  single `DPI` constant, a `theme_set`/`theme_update` design system with `SZ_*`
  size constants, and any analysis/plot helpers folded in so the figure is
  self-contained.
- **`panels/panel<A…>.R`** — each panel sources setup, loads only the inputs it
  needs, builds one `ggplot`/`patchwork` object, and caches it to `objects/`.
- **`assemble.R`** — composes the cached panel objects (via
  [patchwork](https://patchwork.data-imaginist.com/)) and writes the final
  `Figure<N>.pdf` plus a PNG preview. Panel tagging happens here and only here.
- **`run_all.R`** — panels → assemble, end to end.
- **`supp_*.R`** — standalone supplemental figures.

**Manifest-driven inputs** (`common/manifest.R`, used by Fig1 and Fig3). Each
figure's `inputs.tsv` maps a logical `key` to one or more files: a
`local_subpath` (relative to a staging root that defaults to the working
directory) and, where known, a `gs://` origin. A file is fetched with `gsutil`
**only when the local copy is missing**, so pre-staged inputs skip the copy.
Panels reach inputs through `input_path(key)`, `load_input(key)` and
`input_files(key)`, and each panel declares only what it needs. To add or
re-point an input you edit `inputs.tsv` — no script changes. (Fig2 predates this
and reads its inputs directly; it has no `inputs.tsv`.)

**Rendering.** Text uses the **Inter** font, fetched from Google Fonts at runtime
via [showtext](https://cran.r-project.org/package=showtext), with an automatic
fallback to the default sans font when Google Fonts is unreachable. Every output
DPI is matched to the `showtext` DPI so text scales correctly.

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

See each figure's `README.md` for its exact panel set, inputs, outputs, and any
figure-specific steps (Fig1's MetaSoft supplemental, Fig3's derived-data cache,
etc.).

### Tests

`tests/` holds unit and regression checks that need **no controlled-access
data** — they run against mock inputs built by `tests/fixtures.R`:

```bash
Rscript tests/run_all.R
```

The suite covers the shared manifest loader and the Fig3 pipeline, including
regression tests for two defects in the original Fig3 code and a test that
asserts the **rendered** panel-letter sequence of the assembled figure.

---

## Requirements

- **R** ≥ 4.3
- Core CRAN packages used across figures: `tidyverse`, `data.table`, `patchwork`,
  `scales`, `showtext`, `sysfonts`, `paletteer`.
- Additional per figure — e.g. `igraph` + `broom` (Fig3), `ggrepel` +
  `rcartocolor` + `colorspace` + `readxl` (Fig1), and the Bioconductor package
  `rtracklayer` (Fig1). Each figure's `README.md` lists its own dependencies and
  an install snippet.
- **`gsutil`** (Google Cloud SDK) for on-demand input staging from the secure
  bucket — present by default on the *All of Us* Workbench.
- External tools where noted: **MetaSoft** (Java) for Fig1's supplemental;
  **FastENLOC** produces Fig3's inputs upstream (Fig3 treats them as inputs, not
  as something it runs).

---

## Data & inputs

Every input derives from *All of Us* controlled-access data staged from a secure
Google Cloud Storage bucket on the Researcher Workbench. For figures with a
manifest (`Fig1/inputs.tsv`, `Fig3/inputs.tsv`), that file is also the record of
each input's origin; not every input has a `gs://` URI yet, and those are
flagged in the corresponding figure README as needing to be pre-staged.

Access to the underlying data is governed by the *All of Us* Data Use Agreement
and is **not** provided through this repository. The pipelines here have been
exercised against unit-test mock fixtures; the real end-to-end run happens on the
Workbench, where the data resides.

---

## Development notes

Fig1 and Fig3 were refactored from single flat notebook exports into the modular
pipelines above. The design specs and step-by-step implementation plans for those
cleanups live under [`docs/superpowers/`](docs/superpowers/), and each cleaned
figure's `NOTES.md` records deferred work and any deliberate difference from the
original analysis (including, for Fig3, several changes that affect figure
output — read it before comparing against a previously published version).
