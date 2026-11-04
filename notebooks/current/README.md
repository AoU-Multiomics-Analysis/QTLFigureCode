# November 2 figure notebooks

The five main notebooks use the November 2, 2026 submission as their reference.
They rebuild current plots from saved results. They also show the exact submitted
image in a separate reference cell. This image is not a new analysis result.

Figure 1 now uses the post-QC schematic, circos plot, seven gene categories,
three-assay credible-set curves, and combined lead-PIP plot. Figure 2 uses the
compact A-F layout, four frequency bins, and REM labels. Figure 3 includes the
current annotation panel E. Figure 4 uses the Science layout with 68 trait labels.
Figure 5 uses the current five panels. ACVRL1 is now in S31; SPR is in S32.

## Run

Use the R kernel. Start in the repository root. Set these variables before you
run a notebook. Each value is a local directory path:

```r
Sys.setenv(
  QTL_PAPER_ROOT = "/path/to/local/paper",
  QTL_SOURCE_ROOT = "/path/to/figure/checkout"
)
```

`QTL_PAPER_ROOT` contains the local `output/` and `Supplemental Figures/` folders.
`QTL_SOURCE_ROOT` contains staged figure data and saved `FigN/objects/` and
`FigN/data/` caches. It can be the current checkout. Input files are not copied
into Git. The R code checks that each requested file exists and is readable.
New exports use the ignored `notebook_outputs/` folder. Set `QTL_NOTEBOOK_OUTPUT`
to use a different output folder.

Run the cells from top to bottom. The notebooks use tidyverse and patchwork.
Some plots also use ggbeeswarm, ggdist, ggh4x, ggrepel, cowplot, showtext, or svglite.
The R batch executor uses jsonlite and evaluate. The R kernel uses IRdisplay.
The recorded `qtl_execution` metadata lists input file hashes and the R version.
It confirms that the input files were unchanged during the recorded execution.
The [package version table](session_packages.tsv) records the local R environment.

To regenerate notebook code, run:

```sh
python3 tools/build_notebooks.py
python3 tools/build_notebooks.py Fig1 S11
```

Generation retains saved outputs only while every preceding code cell is
unchanged. A changed setup or code cell clears its outputs and all later outputs.
Execute the notebook again after any code change. To clear all outputs, add
`--clear-outputs`.

To execute with the local R environment without a Python Jupyter installation:

```sh
Rscript tools/execute_current_notebooks.R
python3 tools/check_current_notebooks.py --require-executed
```

The original pipeline generator remains available with `--legacy`. The existing
`SuppFig1` notebook remains available for its Workbench inputs. It was not
re-executed because those inputs were not found in this local source collection.

## Supplements

The new notebooks are in `notebooks/supplements/`. They cover S4, S11-S16,
S20-S25, S28, S30, and S31. The [coverage table](coverage.tsv) lists all 32 sections.
It distinguishes complete source data, missing inputs, partial inputs, and the
S27 control-set review. The section number follows the November 2 manuscript;
earlier working drafts used different numbers.

S11-S14 use the current EAS versions. Their population order is Combined, AFR,
AMR, EUR, SAS, EAS, MID. S14 retains the fixed shared-gene set and the saved
descriptive curves. These notebooks do not refit SuSiE models.

S19 has saved event counts but lacks complete model results. S32 has the
fine-mapping panel but lacks the complete Borzoi inputs. S27 has a recorded
difference between the control set in the plot and the caption comparisons.
These three sections do not have new notebooks in this update.

## Scope and checks

The notebooks preserve the recorded scientific results. They do not repeat
association mapping, variant annotation, or model fitting unless the code
section says so. Some sections restore a saved plot object. Others calculate
summaries or rebuild plots from saved tables. Each section states its scope.
Rendering evaluates ggplot statistics, such as densities and straight-line fits.
This does not repeat the upstream association or fine-mapping analyses.
The complex Figure 1 artwork and the final vector label positions are retained
in the exact submission reference. The R reconstructions can have different
spacing or font metrics.

The local check executes each supported notebook, checks source file hashes,
and allows aggregate images and status messages in outputs. The GitHub Actions
smoke test checks the generated code, saved images, and R syntax. It does not
have controlled data and does not rerun the scientific analyses.
