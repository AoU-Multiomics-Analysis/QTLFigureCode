# Figure 2 — causal-variant genetics across expression & splicing QTLs

Code to generate **Figure 2** (and one supplemental figure) for the manuscript.
The figure summarizes properties of SuSiE fine-mapped causal variants across
molecular QTL assays in the *All of Us* cohort: their allele-frequency spectrum,
ancestry enrichment, molecular effect sizes, and their relationship to
expression-outlier and dosage-like phenotypes. Most — but not all — panels
restrict to the PIP > 0.9 subset; see
[Fine-mapping filter](#fine-mapping-filter) below.

> **Data note.** The inputs are derived from *All of Us* **controlled-access**
> genomic data and are **not** included in this repository. See
> [Data & inputs](#data--inputs). Do not commit the staged inputs or the derived
> `.rds`/figure outputs to a public repository.

---

## The figure

**Main figure** (`Figure2.pdf` / `Figure2.png`), composed of seven panels (A–G):

| Panel | Content | Source |
|-------|---------|--------|
| **A** | Cumulative MAF distribution of causal variants (fraction with MAF ≤ k) | expression + splicing (proteomics when available) |
| **B** | Causal-variant population-enrichment distributions by primary population, rare vs. common, with ancestry-enriched variant counts | **all** fine-mapped variants (PIP > 0.9), pooled across assays, one row per variant |
| **C** | \|log2 aFC\| effect size by variant consequence, rare (dark) vs. common (light), annotated with the per-consequence rare−common effect shift | expression |
| **D** | Expected heterozygous expression from cumulatively combining a gene's independent variants (deletion-like / duplication-like) | expression (aFC) — **all** variants with a defined `log2_aFC`, *not* the PIP > 0.9 subset |
| **E** | Expression-outlier enrichment across effect-size bins (underlier vs. overlier) | QTL-burden summary |
| **F** | Median dosage-affected genes per individual by percent-change bin (protein-coding) | QTL-burden summary |
| **G** | Validated dosage-like events (deletion- vs. duplication-like) | QTL-burden counts |

**Supplemental figure** (`FigS_AncestryStacked.pdf` / `.png`): composition of
ancestry-enriched variants across AoU allele frequency (stacked by enriched
ancestry), with the proportion of variants shared across all continents overlaid
as a line. Computed over fine-mapped expression variants (PIP > 0.9).

### Fine-mapping filter

The PIP > 0.9 cutoff is applied **per panel**, not once upstream, and panel D
does not apply it:

- Panels **A**, **B** and **C** filter on `pip > 0.9` in `01_prep_data.R`
  (lines ~29, ~42 and ~79), as does `supp_ancestry_stacked.R` (line ~21).
- Panel **D** starts from the unfiltered `eQTL_COMB_susie` table and filters
  only on `!is.na(log2_aFC)` (`01_prep_data.R` line ~121) — i.e. it is computed
  over *every* expression variant with a defined allelic fold change, not just
  the fine-mapped ones.

*Open question:* whether panel D is meant to use the full expression table or
the PIP > 0.9 subset is worth confirming with the analysis author; the code is
documented here as-is.

---

## Repository layout

```
.
├── common/manifest.R            # shared input manifest + gsutil staging (all figures)
├── fig2/                        # modular pipeline (one directory per figure)
│   ├── 00_setup.R               #   libraries, font, theme, uniform text sizes, manifest loaders, shared helpers
│   ├── inputs.tsv               #   input manifest: key, Assay, local_subpath, gs_uri
│   ├── 01_prep_data.R           #   reads raw inputs ONCE (by manifest key) -> fig2/data/<A-G>.rds
│   ├── panels/panelA.R … panelG.R  #   each builds one panel -> fig2/objects/<X>.rds
│   ├── assemble.R               #   composes panels -> fig2/Figure2.pdf / .png
│   ├── run_all.R                #   prep -> panels -> assemble (end-to-end)
│   ├── supp_ancestry_stacked.R  #   supplemental figure (standalone)
│   ├── NOTES.md                 #   behaviour changes + deferred work
│   ├── Fig2.ipynb               #   generated notebook (tools/build_notebooks.py)
│   ├── data/                    #   per-panel data frames (.rds)   [git-ignored]
│   ├── objects/                 #   per-panel ggplot objects (.rds) [git-ignored]
│   ├── Figure2.pdf / Figure2.png   # main figure outputs           [git-ignored]
│   └── FigS_AncestryStacked.*   #   supplemental figure outputs    [git-ignored]
├── COMB.*_SusieMerged.annotated.ancestry_skew.tsv.gz   # staged SuSiE inputs [not committed]
└── Fig2Data/                    # staged burden inputs (controlled access) [not committed]
```

There is one pipeline, not two. The large SuSiE tables are read **once** in
`01_prep_data.R`; each panel is then built and cached independently, and
`assemble.R` composes them — so you can iterate on a single panel without
re-reading the raw data. Per-panel **data** and **ggplot objects** are both
cached, so a downstream assembly/re-plot script can consume either. Both caches
are participant-derived and git-ignored.

Note that the directory on disk is `Fig2/`, but the scripts self-reference as
lowercase `fig2/` — as everywhere else in this repo.

---

## Running

Run everything from the **project root** (the folder the inputs are staged
under).

```bash
# Full pipeline: data prep -> panels -> assembly
Rscript fig2/run_all.R

# Supplemental figure (standalone)
Rscript fig2/supp_ancestry_stacked.R
```

Rebuild a single panel after editing it (no need to re-read the raw data):

```bash
Rscript fig2/panels/panelC.R   # rebuild panel C from its cached data
Rscript fig2/assemble.R        # recompose the figure
```

Force a full rebuild of the per-panel data and objects:

```bash
rm -rf fig2/data fig2/objects && Rscript fig2/run_all.R
```

### Outputs

Written into `fig2/` by default (`FIG2_OUTDIR`):

- `Figure2.pdf` — 40×24″ vector figure (primary deliverable)
- `Figure2.png` — full-size raster preview
- `FigS_AncestryStacked.pdf` / `.png` — supplemental figure

---

## Requirements

- **R** ≥ 4.3
- CRAN packages: `tidyverse`, `data.table`, `scales`, `patchwork`, `ggbeeswarm`,
  `ggdist`, `ggh4x`, `paletteer`, `rcartocolor`, `colorspace`, `showtext`,
  `sysfonts`, `slider`, `rlang`, `broom`

```r
install.packages(c(
  "tidyverse","data.table","scales","patchwork","ggbeeswarm","ggdist","ggh4x",
  "paletteer","rcartocolor","colorspace","showtext","sysfonts","slider","rlang","broom"
))
```

Text uses the **Inter** font, fetched from Google Fonts at runtime via
`showtext`; if it can't be reached, the scripts fall back to the default sans
font automatically. All output DPI is matched to the `showtext` DPI so text
sizes render correctly.

---

## Data & inputs

**Inputs are manifest-driven.** `fig2/inputs.tsv` maps each logical input (a
`key`) to a `local_subpath` (relative to the staging root, `FIG2_INPUT`, which
defaults to the cwd) and, where applicable, a `gs://` origin. `01_prep_data.R`
reads each input by key via helpers in `00_setup.R`:

- `read_qtl(key)` → `fread(input_path(key))` for a single-file key (staged if missing)
- `input_path(key)` → local path for a single-file key
- `load_input(key)` → reads + row-binds every file for a key, attaching the
  manifest's `Assay` metadata column

When a file is absent locally and the manifest has a `gs://` URI, the loader
copies it down with `gsutil` **once**, into `<cwd>/<local_subpath>` (the local
copy then acts as the cache); pre-staged files skip the copy entirely. All six
Fig2 inputs are locally-derived (`gs_uri = -`), so they must already be present
under the staging root — the three SuSiE tables at its top level and the three
QTL-burden inputs under `Fig2Data/`. There is no remote origin recorded for
them.

| key | File (`local_subpath`) | Assay | Used by |
|-----|------------------------|-------|---------|
| `susie_eqtl` | `COMB.expression_SusieMerged.annotated.ancestry_skew.tsv.gz` | Expression | A, B, C, D |
| `susie_sqtl` | `COMB.splicing_SusieMerged.annotated.ancestry_skew.tsv.gz` | Splicing | A, B |
| `susie_pqtl` | `COMB.proteomics_SusieMerged.annotated.ancestry_skew.tsv.gz` | Proteomics | A, B *(optional — auto-included if present)* |
| `burden_enrichment` | `Fig2Data/QTLBurdenOutlierEnrichment.tsv` | — | E |
| `median_ko` | `Fig2Data/MedianKOSummary.tsv` | — | F |
| `burden_count` | `Fig2Data/QTLBurdenOutlierCount.rds` | — | G |

The expression/splicing tables are SuSiE fine-mapping results annotated with
variant consequences, allele-count/-frequency data (per-ancestry GVS fields),
allelic fold change (`log2_aFC`, expression only), and variant-level
ancestry-skew statistics (`gvs_odds_ratio`, `gvs_p_value`). Proteomics is
optional — the scripts build the figure without it and fold it in automatically
if the file appears.

Access to the underlying *All of Us* data is via the *All of Us* Researcher
Workbench under its Data Use Agreement; it is not shared through this repository.
