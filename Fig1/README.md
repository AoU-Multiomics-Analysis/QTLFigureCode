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
│   ├── 00_setup.R               #   libraries, font, theme, sizes, palette, folded-in functions
│   ├── 01_prep_data.R           #   reads raw inputs ONCE -> fig1/data/<B-H>.rds
│   ├── panels/panelA.R … panelH.R  #   each builds one panel -> fig1/objects/<X>.rds
│   ├── assemble.R               #   composes panels -> Figure1.pdf / Figure1.png
│   ├── run_all.R                #   prep -> panels -> assemble (end-to-end)
│   ├── supp_metasoft.R          #   supplemental MetaSoft figure (standalone)
│   ├── download_inputs.sh       #   provenance record of gs:// inputs (not executed)
│   ├── NOTES.md                 #   deferred analyses (novel pGene) + notes
│   ├── data/                    #   per-panel data frames (.rds)   [git-ignored]
│   └── objects/                 #   per-panel ggplot objects (.rds) [git-ignored]
├── input_data/                  # raw inputs (controlled access, NOT committed)
└── Figure1.pdf / Figure1.png    # main figure outputs
```

The pipeline is cached at each boundary: `01_prep_data.R` reads the (large) raw
inputs **once**; each panel is built and cached independently; then `assemble.R`
composes them. Iterate on a single panel without re-reading the raw data. Both
per-panel **data** and per-panel **ggplot objects** are saved.

The analysis/plot functions were folded in from `QTLAnalysisUtils/QTLGeneAnalysis`
so this figure is self-contained.

---

## Running

Run everything from the **project root** (the folder containing `input_data/`).

```bash
# Modular pipeline (recommended): data prep -> panels -> assembly
Rscript fig1/run_all.R

# Supplemental MetaSoft figure (requires METASOFT_JAR)
METASOFT_JAR=/path/to/Metasoft.jar Rscript fig1/supp_metasoft.R
```

Rebuild a single panel after editing it (no need to re-read the raw data):

```bash
Rscript fig1/panels/panelF.R   # rebuild panel F from its cached data
Rscript fig1/assemble.R        # recompose the figure
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

Expected under `input_data/` (not distributed here). See
[`download_inputs.sh`](download_inputs.sh) for the `gs://` provenance of each file.

| File / group | Used by |
|--------------|---------|
| `reference/gencode.v48.annotation.gtf.gz` | D (gene types) |
| `reference/echo_v4_r2.ancestry_preds.tsv` | B, C |
| `reference/COMB_eQTL_9k_sample_list.tsv` | B |
| `permutations/{eqtl,sqtl,pqtl}/*.cis_qtl.txt.gz` (7 ancestries each) | C, D |
| `finemapping/eqtl/AoU_*_9k_eQTL_SusieMerged.tsv.gz` | E, G |
| `finemapping/sqtl/AoU_*_9k_sQTL_SusieMerged.annotated.tsv` | G |
| `finemapping/pqtl/pQTL_fm_data_10k.tsv` | F, G |
| `annotated/AoU_COMB_9k_eQTL_SusieMerged.annotated.tsv` | H |
| `annotated/AoU_COMB_9k_sQTL_SusieMerged.annotated.tsv` | H |
| `annotated/pQTL_susie_annotated.tsv` | H |
| `external/Whole_Blood.v10.eGenes.txt.gz` (GTEx) | D |
| `external/Whole_Blood.v10.sGenes.txt.gz` (GTEx) | D |
| `external/INTERVAL_eQTL_summary.xlsx` | D |
| `external/INTERVAL_cis_sQTL_summary.xlsx` | D |
| `external/Whole_blood.maf0001.100.permutations.txt` (TOPMed) | D |
| *UKBB pQTL reference* | D (novel pGene — **deferred**, see NOTES.md) |

> **Note.** eQTL and sQTL permutation files share basenames
> (`AFR.cis_qtl.txt.gz`, …); keep them in the separate `permutations/eqtl/` and
> `permutations/sqtl/` subdirectories so they don't collide.

Access to the underlying *All of Us* data is via the *All of Us* Researcher
Workbench under its Data Use Agreement; it is not shared through this repository.
