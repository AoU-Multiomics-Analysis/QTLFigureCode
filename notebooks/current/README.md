# Figure notebooks

The five main notebooks and 16 supplement notebooks read local result tables,
show the data processing, and build their plots in R. Their cells include the
filters, identifier matching, joins, counts, statistical models, and plotting
functions. Intermediate objects remain in the current R session. The notebooks
do not read saved data or plot RDS files.

The scientific reference is the November 2, 2026 figure set. The main notebooks
also display the complete reference PNG after their rebuilt panels.

## Run

Use the R kernel. Start in the repository root. Set the local input directories:

```r
Sys.setenv(
  QTL_INPUT_ROOT = "/path/to/source/result/tables",
  QTL_SUSIE_ROOT = "/path/to/separate/fine-mapping/tables",
  QTL_PAPER_ROOT = "/path/to/paper/references",
  QTL_SOURCE_ROOT = "/path/to/figure/checkout"
)
```

`QTL_INPUT_ROOT` contains the staged result-table paths in each figure's
`inputs.tsv`. `QTL_SUSIE_ROOT` contains the three combined-assay fine-mapping
exports. `QTL_PAPER_ROOT` contains the frozen annotations, external study results,
aggregate source counts, and reference images. `QTL_SOURCE_ROOT` supplies public
reference tracks and source metadata. These roots are read in place.

Run the cells from top to bottom. The notebook checks each input file and records
its hash. New exports use the ignored `notebook_outputs/` directory. Set
`QTL_NOTEBOOK_OUTPUT` to use a different export directory. Only aggregate plots
and status messages are retained in notebook outputs.

The notebooks use tidyverse and patchwork. Individual figures also use packages
such as arrow, data.table, ggbeeswarm, ggdist, ggh4x, ggrepel, cowplot, rcartocolor,
MatchIt, slider, and jsonlite. The [package table](session_packages.tsv) records the local
R environment. The batch executor uses evaluate; the R kernel uses IRdisplay.

## Input versions

The [preparation guide](../../preparation/README.md) and
[manifest](../../preparation/manifest.json) list source paths and file overrides.
For S11-S14, `QTL_FIG1_EAS_EQTL_FILE` must identify the current EAS eQTL result.
The older similarly named file is a different analysis. Its current download is
gzip CSV despite its `.tsv` suffix; the reader checks the file signature.

Use `QTL_TOPMED_CS_FILE` for the TOPMed credible-set result. Use
`QTL_FIG1_PROTEIN_BED_FILE` and `QTL_FIG1_PROTEIN_COVARIATES_FILE` for the protein
sample headers. S24/S25 use `QTL_FIG4_MEMBERSHIP` and `QTL_FIG4_SIGNAL` when the
staged input directory has an older colocalization export.

The notebooks calculate their own intermediate values. `QTL_PREP_ROOT` is used
only by the separate preparation runner; it cannot replace notebook inputs.

## Code and generation

Authored R files are in this directory. A `# @include repository/path.R` line
inserts that file's actual R code into the notebook. Helpers are visible in the
cells. An included-code change invalidates saved outputs from that point onward.
No runtime `source()` call hides a processing step.

```sh
python3 tools/build_notebooks.py
python3 tools/build_notebooks.py Fig1 S11
Rscript tools/execute_current_notebooks.R
python3 tools/check_current_notebooks.py --require-executed
```

Generation retains outputs only while all preceding code cells are unchanged.
Execute again after a code change. Add `--clear-outputs` to remove every output.
The old modular-pipeline generator remains available with `--legacy`.

## Supplements and input gaps

The supplement notebooks cover S4, S11-S16, S20-S25, S28, S30, and S31. The
[coverage table](coverage.tsv) lists all 32 manuscript sections.

Figure 1 pre-QC availability uses supplied aggregate counts because the
participant-availability export is unavailable. Post-QC counts are recalculated
from sample lists and protein headers.

Figure 4A shows its source-processing and plotting code. If the full harmonized
gene-level table is absent, it reports a skip and leaves A empty. B-D still run.
It does not substitute a saved plot.

S30 shows the recovered aggregate count-to-percentile calculation. The original
genome-wide ranking inputs and procedure are unavailable. The supplied formulas
remain its input boundary.

S19 lacks complete model results. S32 lacks complete Borzoi inputs. S27 has a
recorded control-set difference. No new notebooks were added for these sections.
The existing `SuppFig1` notebook remains available for its Workbench inputs.

## Checks

Local verification compares fresh counts, data frames, model results, and plot
layers with the accepted outputs. Input hashes are checked after execution.
GitHub Actions checks generation, visible code, saved aggregate images, R
syntax, and synthetic processing/plot tests without controlled inputs.
