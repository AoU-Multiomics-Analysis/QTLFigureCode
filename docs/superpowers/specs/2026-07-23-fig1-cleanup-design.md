# Figure 1 Cleanup — Design

**Date:** 2026-07-23
**Author:** Evin Padhi (with Claude)
**Status:** Proposed

## Goal

Refactor `Fig1/` from a single 548-line linearized-notebook script (`Fig1.R`) into the
same clean, cached, multi-stage pipeline used by `Fig2/`. Fold the required functions from
`/Users/epadhi/QTLAnalysisUtils` into the directory so Fig1 is self-contained, fix the
latent bugs surfaced during reconciliation, and document the (pre-staged) inputs.

## Source of truth

The authoritative source is the notebook `QTLAnalysisUtils/notebooks/QTL_paper_Fig1.ipynb`,
which `Fig1.R` was exported from. The notebook's per-panel cells are the correct logic. Its
final assembly cell contains the same variable-name bugs as `Fig1.R` (the plot functions'
internal return-variable names were pasted in by mistake); these are resolved as:

- `plot_final` → `expression_plot`
- `novel_xGenes` → `novel_qtl_count_plot`
- `vertical_enrichment_plot` → `EnrichmentPlot`

## Decisions (locked)

| Decision | Choice |
|---|---|
| Source of truth | Reconcile with `QTL_paper_Fig1.ipynb` |
| Structure | Match Fig2 exactly |
| Utils placement | Fold used functions into `00_setup.R` |
| Data access | Assume pre-staged inputs in `input_data/`; no runtime `gsutil`/Java in the main pipeline |
| MetaSoft | Keep as a standalone `supp_metasoft.R` supplemental (not part of Figure 1) |
| UKBB pQTL / novel pGene | Drop the pGene bar from Panel D for now; record the deferred analysis |

## Target directory structure

```
Fig1/
├── README.md               # same section skeleton as Fig2/README.md
├── 00_setup.R              # libs, env-overridable paths, theme, palette, folded-in functions (fixed)
├── 01_prep_data.R          # reads input_data/ once → data/<A–H>.rds bundles
├── run_all.R               # prep → panels → assemble
├── assemble.R              # patchwork compose → Figure1.pdf / Figure1.png
├── supp_metasoft.R         # standalone MetaSoft supplemental → xGeneMetasoftSummary.pdf
├── download_inputs.sh      # commented, NON-executed record of gs:// provenance (controlled-access)
├── NOTES.md                # deferred analyses (e.g. UKBB pGene novelty)
├── panels/
│   ├── panelA.R  # schematic placeholder (blank_plot slot; study-design diagram added externally)
│   ├── panelB.R  # ancestry counts in cohort
│   ├── panelC.R  # proportion of genes with significant QTLs (expression | splicing + protein)
│   ├── panelD.R  # novel eGene + sGene counts vs external consortia (pGene deferred)
│   ├── panelE.R  # CDF of credible sets per eGene (COMB)
│   ├── panelF.R  # credible-set sizes across pQTL ancestry groups (+ UKBB reference point)
│   ├── panelG.R  # proportion of credible sets with a single causal variant
│   └── panelH.R  # fine-mapping functional enrichment by assay
├── data/                   # per-panel data bundles (.rds), git-ignored
├── objects/                # per-panel ggplot objects (.rds), git-ignored
└── input_data/             # pre-staged raw inputs, git-ignored (controlled-access)
    ├── permutations/{eqtl,sqtl,pqtl}/   # subdir'd to avoid basename collisions (see note)
    ├── finemapping/{eqtl,sqtl,pqtl}/
    ├── annotated/                        # COMB annotated susie tables for enrichment
    ├── external/                         # GTEx / INTERVAL / TOPMed reference tables
    └── reference/                        # gencode gtf, ancestry preds, sample list
```

Panel letters follow the notebook's `plot_annotation(tag_levels = 'A')` order. Panel A is a
blank slot (as in Fig2's `blank_plot`) where the externally-drawn study-design schematic is placed.

## Pipeline stages (mirrors Fig2)

### `00_setup.R` — shared foundation (cheap, idempotent)
- All `library()` calls in `suppressPackageStartupMessages({...})`:
  `tidyverse, data.table, rtracklayer, plyranges, rcartocolor, scales, ggrepel, patchwork,
  readxl, showtext, sysfonts, paletteer, colorspace`.
- Path constants overridable via env vars, e.g.
  `INPUT <- Sys.getenv("FIG1_INPUT", "input_data")`, plus `DATA_DIR`, `OBJ_DIR`, `OUTDIR`;
  `dir.create()` the output dirs.
- Font/DPI: `font_add_google("Inter")` with `"sans"` fallback; `DPI <- 300`;
  `showtext_opts(dpi = DPI)` (comment noting ggsave DPI must match).
- Theme: `theme_set(theme_classic())` + `theme_update(...)`; `AncestryColorPalette`; `blank_plot`.
- **Folded-in functions (used-only, fixed):**
  - Loaders: `LoadAncestryCounts`, `LoadGTExeQTL`, `LoadGTExsQTL`, `LoadIntervaleQTL`,
    `LoadIntervalsQTL`, `LoadTopmedeQTL`, `LoadExternaleQTLs`, `LoadExternalsQTLs`.
  - Transforms: `ExtractNovelQTLs`, `SummarizeSigQTLs`, `ComputeCredibleSetCDF`,
    `ComputeFinemappingAnnotationEnrichment`.
  - Plot builders: `PlotAncestryCounts`, `PlotPropQTLSig`, `PlotNovelQTLCount`, `PlotCSSizes`,
    `PlotFinemappingEnrichment`.
  - Dropped as unused by Figure 1: `FormatMetaSoft`, `PlotMetaSoftMvaluesProps`,
    `PlotMetaSoftMvalueByPop`, `LoadUKBBpQTL` (MetaSoft trio moves to `supp_metasoft.R`).

### `01_prep_data.R` — all wrangling, one bundle per panel
Reads `input_data/` once, computes shared intermediates (`QTL_data` from the permutations,
`gene_types` from the GTF, `ancestry_counts`, the fine-mapping tables), then writes a named-list
`.rds` per panel via a `saveit(x, name)` helper that logs its path.

| Bundle | Contents | Raw inputs |
|---|---|---|
| `data/B.rds` | `ancestry_counts` | ancestry preds + COMB sample list |
| `data/C.rds` | `expression_only`, `splicing_proteomics` (via `SummarizeSigQTLs`) | all QTL permutations, ancestry_counts, palette |
| `data/D.rds` | `novel_qtl_count` (eGene + sGene) | QTL_data COMB, gene_types, `ExternaleQTLs`, `ExternalsQTLs` |
| `data/E.rds` | `ccdf` (via `ComputeCredibleSetCDF`) | eQTL fine-mapping |
| `data/F.rds` | `pQTL_cs_sizes`, `group_means`, `common_var_means`, `ref_means` | pQTL fine-mapping (`pQTL_fm_data_10k.tsv`) |
| `data/G.rds` | `summary_data`, `causal_variant_count` | eQTL/sQTL/pQTL fine-mapping |
| `data/H.rds` | `QTLFinemappingEnrichments`, `variant_functional_groups` | annotated susie tables (e/s/pQTL COMB) |

(Panel A needs no data bundle.)

### `panels/panelX.R` — thin, plot-only
Read one `data/X.rds`, call the relevant `Plot*` builder from setup (or a small inline ggplot for
E/G which are bespoke), apply panel-specific theme tweaks, `saveRDS()` to `objects/X.rds`, and end
with a `message("Panel X built.")`.

### `assemble.R`
Load `objects/*.rds` into a named list, reproduce the `top_row / second_row / third_row / bottom_row`
patchwork layout from the notebook, `plot_annotation(tag_levels = 'A')`, save `Figure1.pdf` (25×27")
then a PNG preview at `showtext_opts(dpi = 100)` (restoring DPI after), matching Fig2's raster pattern.

### `run_all.R`
`source("00_setup.R")` → `source("01_prep_data.R")` → loop `panels/panel{A..H}.R` → `source("assemble.R")`.

### `supp_metasoft.R` (supplemental)
Standalone script keeping `FormatMetaSoft` + the two MetaSoft plot builders (fixed), producing
`xGeneMetasoftSummary.pdf`. Documented external dependency: the MetaSoft Java tool and its
`.metasoft.tsv` output (either pre-staged or produced by a documented manual step). Not part of the
main `run_all.R`.

## Bug fixes applied during the fold-in

1. `LoadExternalsQTLs`: bind `GTEx_sQTL_df` (was the undefined `GTEx_eQTL_df`).
2. Remove the duplicate `LoadGTExeQTL` definition.
3. `SummarizeSigQTLs`: join on its `AncestryCounts` argument, not the global `ancestry_counts`.
4. `ComputeFinemappingAnnotationEnrichment` / `PlotFinemappingEnrichment`: take their inputs as
   arguments instead of reading module globals (`sQTL_/eQTL_/pQTL_variant_enrichment`,
   `variant_functional_groups`).
5. Assembly variable-name fixes (`plot_final`/`novel_xGenes`/`vertical_enrichment_plot`).
6. `ptions(...)` → `options(...)` typo (Fig1.R only).
7. `casual` → `causal` throughout Panel G (`causal_variant_count`, `number_causal`, labels).
8. Permutation basename collision: eQTL and sQTL permutation files share basenames
   (`{ANCESTRY}.cis_qtl.txt.gz`); the old code downloaded both to the same dir. Pre-staged inputs
   are subdir'd (`permutations/eqtl/`, `permutations/sqtl/`, `permutations/pqtl/`) to keep them distinct.

## Inputs (pre-staged; documented in README "Data & inputs")

- **Permutations** (7 ancestries × 3 assays): tensorQTL `cis_qtl.txt.gz` per ancestry.
- **Fine-mapping**: eQTL & sQTL SusieMerged per ancestry; pQTL `pQTL_fm_data_10k.tsv`.
- **Annotated susie** (COMB): eQTL, sQTL (`*.annotated.tsv`), pQTL `pQTL_susie_annotated.tsv`.
- **External references**: GTEx v10 Whole Blood e/sGenes; INTERVAL e/sQTL xlsx; TOPMed Whole Blood eQTL.
- **Reference**: `gencode.v48.annotation.gtf.gz`, `echo_v4_r2.ancestry_preds.tsv`, COMB sample list.
- **Deferred** (see NOTES.md): UKBB pQTL supplementary table for the novel-pGene comparison.

The exact `gs://` origins of every file are recorded, commented and non-executed, in
`download_inputs.sh` so provenance is preserved without coupling the pipeline to the Terra/GCS
environment.

## Deferred / recorded (NOTES.md)

- **Novel pGene novelty (Panel D):** requires a UKBB pQTL reference (`UKBB_pqtl_path`, currently
  unassigned upstream). Once supplied via `input_data/external/`, re-enable `novel_pQTLs` in
  `01_prep_data.R` and add the pGene bar back to Panel D. The plotting path already supports it.

## Out of scope

- No changes to Fig2 or other figure directories.
- No re-running of the analysis against live GCS data; the pipeline consumes pre-staged inputs.
- No new scientific analysis beyond faithfully reproducing the notebook's Figure 1.
