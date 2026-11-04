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
| **D** | Novel eGene / sGene / pGene counts vs. external consortia, by gene type | QTL permutations (COMB) + GTEx / INTERVAL / TOPMed (e/sGenes) + UKBB-PPP (pGenes) |
| **E** | CDF of credible sets per eGene (≥ k credible sets) | eQTL SuSiE fine-mapping (COMB) |
| **F** | Credible-set sizes across pQTL ancestry groups, with a UKBB reference point | pQTL SuSiE fine-mapping |
| **G** | Proportion of credible sets with a single causal variant, by assay | e/s/pQTL SuSiE fine-mapping (COMB) |
| **H** | Fine-mapping functional enrichment of causal variants, by assay | annotated SuSiE tables (COMB) |

The novel-**pGene** bar in Panel D is compared against the UK Biobank Pharma
Proteomics Project (UKBB-PPP) cis-pQTL reference; lncRNA is excluded from that
comparison because the proteomic assay cannot measure it (`panels/panelD.R`).

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
│   ├── enrichment_glm.R         #   alternative enrichment: glm(annotated ~ PIP) (standalone)
│   ├── enrichment_glm_plot.R    #     its figure; sourced by the above, or run alone
│   ├── NOTES.md                 #   result-changing notes (novel pGene, panel H)
│   └── objects/                 #   per-panel ggplot objects (.rds) [git-ignored]
├── Fig1Data/                    # most inputs stage here under the cwd (controlled access, NOT committed)
│   ├── eQTL/ pQTL/ sQTL/        #   permutation results, one file per ancestry
│   ├── PopulationLevelFinemapping/{eQTL,pQTL,sQTL}/   #   SuSiE fine-mapping
│   └── …                        #   external references (GTEx / INTERVAL / TOPMed / UKBB-PPP)
├── reference/                   # COMB_eQTL_9k_sample_list.tsv only
├── gencode.v48.annotation.gtf.gz          # staged at the root, not under Fig1Data/
├── echo_v4_r2.ancestry_preds.tsv          #   "
├── COMB.*_SusieMerged.annotated.ancestry_skew.tsv.gz  # panel H; shared with Fig2/Fig3
└── Figure1.pdf / Figure1.png    # main figure outputs
```

**Inputs are manifest-driven.** `fig1/inputs.tsv` maps each logical input (a
`key`) to one or more files — a `local_subpath` (relative to the staging root,
which defaults to the current working directory). The public manifest omits
private GCS addresses. Each panel declares and loads **only the inputs it needs** via
helpers in `00_setup.R`:

- `input_path(key)` → local path for a single-file key (staged if missing)
- `load_input(key)` → reads + row-binds every file for a key, attaching
  `group` / `Assay` from the manifest
- `input_files(key)` → the manifest rows plus resolved local `path`s

Supply the private input files at their listed local paths before running.
A separate local manifest can retain the original addresses; see the
[private manifest guide](../docs/private-input-manifests.md).
When a file is absent locally and that manifest has a `gs://` URI, the loader
copies it down with `gsutil` **once**, into `<cwd>/<local_subpath>` (the local
copy then acts as the cache); pre-staged files skip the copy entirely. Staged
inputs therefore appear under the working directory, mostly beneath `Fig1Data/`
— keeping their subpaths (rather than a flat dump) so the eQTL and sQTL
permutation files, which share basenames, don't collide. A few keys resolve to
bare filenames at the staging root instead (the GENCODE GTF, the ancestry
predictions, and the three annotated SuSiE tables). The
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

# Alternative annotation enrichment: logistic regression on PIP, standalone
# figure. Needs the three annot tables; takes a few minutes on the eQTL table.
Rscript fig1/enrichment_glm.R
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
- `Figure1_enrichment_glm.pdf` / `.png` — standalone PIP-regression enrichment,
  with its fitted slopes in `fig1/objects/enrichment_glm.tsv`

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

Every input is declared in **`fig1/inputs.tsv`** (the manifest).
Columns: `key`, `group`, `Assay`,
`local_subpath` (relative to the staging root — the cwd by default, overridable
via `FIG1_INPUT`), `gs_uri` (`-` = no public download address).
Supply the files at `<FIG1_INPUT>/<local_subpath>`, or use the working directory
when `FIG1_INPUT` is unset.

| Manifest key | Files | Used by |
|--------------|-------|---------|
| `gene_gtf` | `gencode.v48.annotation.gtf.gz` (staging root) | D (gene types) |
| `ancestry_preds` | `echo_v4_r2.ancestry_preds.tsv` (staging root) | B, C |
| `comb_sample_list` | `reference/COMB_eQTL_9k_sample_list.tsv` | B |
| `perm` | `Fig1Data/{eQTL,sQTL,pQTL}/*.cis_qtl.txt.gz` (7 ancestries × 3 assays) | C, D |
| `finemap_eqtl` | `Fig1Data/PopulationLevelFinemapping/eQTL/AoU_*_9k_eQTL_SusieMerged.tsv.gz` (7) | E, G |
| `finemap_sqtl` | `Fig1Data/PopulationLevelFinemapping/sQTL/AoU_*_9k_sQTL_SusieMerged.annotated.tsv` (7) | G |
| `finemap_pqtl` | `Fig1Data/PopulationLevelFinemapping/pQTL/pQTL_fm_data_10k.tsv` | F, G |
| `annot` | `COMB.{expression,splicing,proteomics}_SusieMerged.annotated.ancestry_skew.tsv.gz` (staging root; shared with Fig2/Fig3) | H |
| `gtex_egenes` / `gtex_sgenes` | `Fig1Data/Whole_Blood.v10.{e,s}Genes.txt.gz` (GTEx) | D |
| `interval_eqtl` / `interval_sqtl` | `Fig1Data/INTERVAL_{eQTL,cis_sQTL}_summary.xlsx` | D |
| `topmed_eqtl` | `Fig1Data/Whole_blood.maf0001.100.permutations.txt` (TOPMed) | D |
| `ukbb_pqtl` | `Fig1Data/41586_2023_6592_MOESM3_ESM.xlsx` (UKBB-PPP, Nature 2023 supplement) | D (novel pGene) |

> **Note.** eQTL and sQTL permutation files share basenames
> (`AFR.cis_qtl.txt.gz`, …); the manifest keeps them in separate
> `Fig1Data/eQTL/` and `Fig1Data/sQTL/` subdirectories so they don't collide.

> **Note.** The three `annot` tables are the same SuSiE exports Figures 2 and 3
> read, so all three figures share one fine-mapping provenance chain; they are
> not Fig1-specific files.

To add or re-point an input, edit `fig1/inputs.tsv` — no script changes needed.

Access to the underlying *All of Us* data is via the *All of Us* Researcher
Workbench under its Data Use Agreement; it is not shared through this repository.
