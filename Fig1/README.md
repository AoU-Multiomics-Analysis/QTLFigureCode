# Figure 1 — multi-omic QTL discovery across ancestries in *All of Us*

Code to generate **Figure 1** (and one supplemental figure) for the manuscript.
The figure summarizes cis-QTL discovery across three molecular assays —
expression (eQTL), protein (pQTL), and splicing (sQTL) — in the ancestrally
diverse *All of Us* cohort: cohort composition, the proportion of genes with
significant QTLs by sample size, novelty relative to prior consortia, and the
fine-mapping (SuSiE) architecture of causal variants.

> **Data note.** The inputs are derived from *All of Us* **controlled-access**
> genomic data and are **not** included in this repository. See
> [Data & inputs](#data--inputs). Do not commit `input_data/` or the derived
> `.rds`/figure outputs to a public repository.

---

## The figure

**Main figure** (`Figure1.pdf` / `Figure1.png`), composed of eight panels:

| Panel | Content | Source |
|-------|---------|--------|
| **A** | Study-design schematic (placeholder slot; diagram added externally) | — |
| **B** | Participant counts per genetic-ancestry group | ancestry predictions + COMB sample list |
| **C** | Proportion of genes with a significant QTL vs. sample size (expression \| splicing + protein) | QTL permutations (all assays × ancestries) |
| **D** | Novel eGene / sGene counts vs. external consortia, by gene type | QTL permutations (COMB) + GTEx / INTERVAL / TOPMed |
| **E** | CDF of credible sets per eGene (≥ k credible sets) | eQTL SuSiE fine-mapping (COMB) |
| **F** | Credible-set sizes across pQTL ancestry groups, with a UKBB reference point | pQTL SuSiE fine-mapping |
| **G** | Proportion of credible sets with a single causal variant, by assay | e/s/pQTL SuSiE fine-mapping (COMB) |
| **H** | Fine-mapping functional enrichment of causal variants, by assay | annotated SuSiE tables (COMB) |

The novel-**pGene** bar in Panel D is deferred pending a UKBB pQTL reference — see
[`NOTES.md`](NOTES.md).

**Supplemental figure** (`xGeneMetasoftSummary.pdf`): cross-ancestry meta-analysis
(MetaSoft) M-value summary — the proportion of genes by number of ancestry groups
supporting the effect, and the M-value distribution per population by assay.
Standalone; requires the MetaSoft Java tool.

---

## Repository layout

```
.
├── fig1/                        # modular pipeline (split by panel)
│   ├── 00_setup.R               #   libraries, font, theme, palette, manifest loaders, folded-in functions
│   ├── inputs.tsv               #   input manifest: key, group, Assay, local_subpath, gs_uri
│   ├── panels/panelA.R … panelH.R  #   each STAGES + LOADS only its own inputs -> fig1/objects/<X>.rds
│   ├── assemble.R               #   composes panels -> Figure1.pdf / Figure1.png
│   ├── run_all.R                #   panels -> assemble (end-to-end)
│   ├── supp_metasoft.R          #   supplemental MetaSoft figure (standalone)
│   ├── NOTES.md                 #   deferred analyses (novel pGene) + notes
│   └── objects/                 #   per-panel ggplot objects (.rds) [git-ignored]
├── permutations/ finemapping/ … # inputs staged here under the cwd (controlled access, NOT committed)
└── Figure1.pdf / Figure1.png    # main figure outputs
```

**Inputs are manifest-driven.** `fig1/inputs.tsv` maps each logical input (a
`key`) to one or more files — a `local_subpath` (relative to the staging root,
which defaults to the current working directory) and, where applicable, a
`gs://` origin. Each panel declares and loads **only the inputs it needs** via
helpers in `00_setup.R`:

- `input_path(key)` → local path for a single-file key (staged if missing)
- `load_input(key)` → reads + row-binds every file for a key, attaching
  `group` / `Assay` from the manifest
- `input_files(key)` → the manifest rows plus resolved local `path`s

When a file is absent locally and the manifest has a `gs://` URI, the loader
copies it down with `gsutil` **once**, into `<cwd>/<local_subpath>` (the local
copy then acts as the cache); pre-staged files skip the copy entirely. Staged
inputs therefore appear under the working directory in `permutations/`,
`finemapping/`, etc. — keeping their subpaths (rather than a flat dump) so the
eQTL and sQTL permutation files, which share basenames, don't collide. The
per-panel **ggplot objects** are cached in `objects/`, so `assemble.R` never
rebuilds them — rebuild a single panel and re-assemble without touching the others.

The analysis/plot functions were folded in from `QTLAnalysisUtils/QTLGeneAnalysis`
so this figure is self-contained.

---

## Running

Run everything from the **project root** (the folder containing `input_data/`).

```bash
# Full pipeline: build every panel (staging inputs as needed) -> assemble
Rscript fig1/run_all.R

# Supplemental MetaSoft figure (requires METASOFT_JAR)
METASOFT_JAR=/path/to/Metasoft.jar Rscript fig1/supp_metasoft.R
```

Rebuild a single panel after editing it (it re-stages/reads only its own inputs):

```bash
Rscript fig1/panels/panelF.R   # rebuild panel F (loads just the pQTL fine-mapping)
Rscript fig1/assemble.R        # recompose the figure from cached objects
```

### Outputs

- `Figure1.pdf` — 25×27″ vector figure (primary deliverable)
- `Figure1.png` — full-size raster preview
- `xGeneMetasoftSummary.pdf` — supplemental figure

---

## Requirements

- **R** ≥ 4.3
- Bioconductor: `rtracklayer` (GTF import)
- CRAN packages: `tidyverse`, `data.table`, `scales`, `patchwork`, `ggrepel`,
  `readxl`, `paletteer`, `rcartocolor`, `colorspace`, `showtext`, `sysfonts`

```r
install.packages(c(
  "tidyverse","data.table","scales","patchwork","ggrepel","readxl",
  "paletteer","rcartocolor","colorspace","showtext","sysfonts"
))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("rtracklayer")
```

Text uses the **Inter** font, fetched from Google Fonts at runtime via
`showtext`; if it can't be reached, the scripts fall back to the default sans
font automatically. All output DPI is matched to the `showtext` DPI so text
sizes render correctly.

External tool (supplemental only): **MetaSoft**
(`http://genetics.cs.ucla.edu/meta/repository/2.0.1/Metasoft.zip`), Java required.

---

## Data & inputs

Every input is declared in **`fig1/inputs.tsv`** (the manifest), which is also the
single record of each file's `gs://` origin. Columns: `key`, `group`, `Assay`,
`local_subpath` (relative to the staging root — the cwd by default, overridable
via `FIG1_INPUT`), `gs_uri` (`-` = local-only). Missing files with a `gs_uri`
are staged under the cwd on demand; pre-stage them at `<cwd>/<local_subpath>`
(or point `FIG1_INPUT` at an existing folder) to skip `gsutil`.

| Manifest key | Files | Used by |
|--------------|-------|---------|
| `gene_gtf` | `reference/gencode.v48.annotation.gtf.gz` | D (gene types) |
| `ancestry_preds` | `reference/echo_v4_r2.ancestry_preds.tsv` | B, C |
| `comb_sample_list` | `reference/COMB_eQTL_9k_sample_list.tsv` | B |
| `perm` | `permutations/{eqtl,sqtl,pqtl}/*.cis_qtl.txt.gz` (7 ancestries × 3 assays) | C, D |
| `finemap_eqtl` | `finemapping/eqtl/AoU_*_9k_eQTL_SusieMerged.tsv.gz` (7) | E, G |
| `finemap_sqtl` | `finemapping/sqtl/AoU_*_9k_sQTL_SusieMerged.annotated.tsv` (7) | G |
| `finemap_pqtl` | `finemapping/pqtl/pQTL_fm_data_10k.tsv` | F, G |
| `annot` | `annotated/AoU_COMB_9k_{e,s}QTL_SusieMerged.annotated.tsv`, `annotated/pQTL_susie_annotated.tsv` | H |
| `gtex_egenes` / `gtex_sgenes` | `external/Whole_Blood.v10.{e,s}Genes.txt.gz` (GTEx) | D |
| `interval_eqtl` / `interval_sqtl` | `external/INTERVAL_{eQTL,cis_sQTL}_summary.xlsx` | D |
| `topmed_eqtl` | `external/Whole_blood.maf0001.100.permutations.txt` (TOPMed) | D |
| *(deferred)* | UKBB pQTL reference | D (novel pGene — see NOTES.md) |

> **Note.** eQTL and sQTL permutation files share basenames
> (`AFR.cis_qtl.txt.gz`, …); the manifest keeps them in separate
> `permutations/eqtl/` and `permutations/sqtl/` subdirectories so they don't collide.

To add or re-point an input, edit `fig1/inputs.tsv` — no script changes needed.

Access to the underlying *All of Us* data is via the *All of Us* Researcher
Workbench under its Data Use Agreement; it is not shared through this repository.
