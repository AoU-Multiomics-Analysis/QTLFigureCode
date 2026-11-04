# Prepare figure data

These scripts produce the inputs used by the plotting notebooks. They include
the source-result filters, identifier matching, annotation joins, counts,
statistical summaries, and the connections to the plot inputs. Input and
generated data files remain outside Git.

The starting points are QTL, fine-mapping, colocalization, TWAS, and reference
result tables. The scripts run locally. They do not submit sequencing,
association-mapping, or cloud fine-mapping jobs.

## Run

Run from the repository root. List the preparation steps and their dependencies:

```sh
python3 tools/prepare_figure_data.py --list
python3 tools/prepare_figure_data.py fig1_novelty --dry-run
```

Set the local input directories with named options:

```sh
python3 tools/prepare_figure_data.py fig1_novelty \
  --input-root /path/to/staged/results \
  --paper-root /path/to/paper \
  --source-root /path/to/figure/checkout \
  --output-root /path/to/prepared
```

`input` contains the original result and reference files. Its paths follow the
figure manifests, such as `Fig1Data/` and the three combined SuSiE tables.
`paper` contains frozen annotations, reference metadata, and accepted outputs
used for comparisons and layout. `source` contains original figure caches and
reference tracks. `susie` can point to a separate source-result directory with
`--susie-root`. The directory names can contain spaces when quoted.

Add `--check-inputs` to check readability before computation. Use a figure name,
such as `Fig1`, to select its steps, or name individual step keys. Dependencies
run first. An unavailable required input stops the run before computation.
Steps with a missing input can be inspected with `--dry-run`.

The output root defaults to `notebook_outputs/prepared/`, which is ignored by
Git. Every generated file is placed in its `paper/` or `source/` namespace.
The runner checks input hashes after each successful step and writes a record
to `logs/`. Outputs cannot replace inputs read by the same step.

## Use the notebooks

The notebooks now perform these calculations in their own cells and construct
their plots from the resulting in-memory tables. They do not read the RDS files
produced by this separate preparation runner. See the
[notebook guide](../notebooks/current/README.md) for their local input roots.

The original runner remains available for standalone preparation and archived
comparison checks. Its connector steps can retain accepted plot templates;
those connectors are not used by the notebooks. The new explicit plot builders
are in `notebooks/current/plots/` and are expanded into the notebook cells.

The [manifest](manifest.json) records the runner's inputs, outputs, dependencies,
original source paths, SHA-256 values, and local checks.

## Input versions and overrides

The EAS eQTL credible-set file must be the source used in the latest analysis.
An older similarly named file has different contents. Supply the exact file
with `--eas-eqtl-file`. The reader handles the gzip CSV export even when its
filename ends in `.tsv`. Use `--topmed-cs-file` for the TOPMed credible-set file.

Use `--protein-bed-file` and `--protein-covariates-file` for the protein cohort
headers when they are outside the staged input tree. S24/S25 use the exact
current colocalization export. Set `--fig4-membership-file` and
`--fig4-signal-file` when the input directory contains an older export.
Equivalent `QTL_*` file variables are listed in the manifest.

## Figure methods

- [Figure 1](Fig1/README.md): source q-value novelty calls, corrected UKBB
  target matching, gene categories, circos coordinates, cohort counts,
  credible sets, PIP summaries, and S4/S11-S13.
- [Figures 2-4](Fig234.md): frequency bins, ancestry tests, aFC and outlier
  summaries, cross-QTL components and models, TWAS summaries, S16/S20-S25,
  and connections to the plot objects.
- [Figure 5](Fig5/README.md): disease-gene and variant selection, ClinVar
  classes, population assignments, the funnel, effect comparisons, Fisher
  tests, S28/S30/S31, and locus plot inputs.
- [S14](S14/prepare_discovery.R): source credible-set counts for seven AoU
  analyses and the external studies, the shared coding-gene intersection,
  participant sample counts, descriptive curves, and GTEx tissue envelopes.
- [S15](S15/prepare_enrichment.R): the original annotation-on-PIP model
  functions in `Fig1/enrichment_glm_fit.R`, P-value adjustment, and plot input.

## Remaining source boundaries

Figure 1 pre-QC availability counts are supplied aggregate inputs. The original
participant-availability export is unavailable. Post-QC sample counts are
recounted from the sample lists and protein file headers.

Figure 4A preparation code is included. Its full harmonized gene-level result
table is unavailable locally, so that scientific step has not run. B-D and the
S24/S25 data were regenerated and compared with the accepted outputs.

S30 retains the recorded aggregate count-to-percentile formulas. Its original
genome-wide ranking procedure and source records are unavailable. The four
display rows can be reproduced from the supplied aggregate inputs.

Some final plot objects use accepted layout or annotation templates. Connector
steps replace their scientific data with prepared values and compare them
before writing the notebook inputs. These dependencies are explicit in the
manifest. The code and the input metadata are uploaded; the templates and
source tables remain local.

## Checks

GitHub Actions checks the manifest, R syntax, notebook generation, and synthetic
filter/count fixtures. It does not have controlled inputs. Local execution
compares regenerated scientific values with accepted outputs. The manifest
records the result of each comparison.
