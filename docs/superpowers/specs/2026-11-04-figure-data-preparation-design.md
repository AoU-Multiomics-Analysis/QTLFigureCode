# Figure data preparation scripts

The user requests the scripts that turn source analysis results into plotted
data, in addition to the plotting notebooks. Source data must remain outside
Git. This extends the approved November 2 figure update in PR 23.

Add preparation code for filtering, identifier matching, annotation joins,
deduplication, counts, derived summaries, and recorded statistical comparisons.
Use the actual local source scripts and the existing modular figure pipelines.
Preserve their scientific rules and document any missing source step. The
boundary is the source analysis-result tables, not raw sequencing or new cloud
association/fine-mapping jobs.

Keep scripts under `preparation/FigN/`. Put shared local-path and output helpers
in `preparation/common.R`. Provide a manifest that maps figure panels to scripts,
required inputs, generated files, dependencies, and original source identities.
A Python runner must accept named local input roots and run scripts in order.
No input table or derived data file is committed. New outputs go below the
ignored `notebook_outputs/prepared/` directory by default.

The preparation outputs use separate `paper/` and `source/` namespaces. The
notebooks can read these through an explicit `QTL_PREP_ROOT` overlay, then fall
back to the original read-only roots for other inputs and reference images.
The current plots must remain the same when no overlay is selected.

Figure 1 must include the complete current novelty filter, corrected UKBB target
mapping, same-assay reference sets, biotype rules, and unique-gene counts. The
other figures and supplements must identify the preparation code for their
current inputs. Uploading scripts that only read the final counts does not meet
this requirement. Existing upstream analysis code can be referenced when it
already supplies the required transformation.

Use tidyverse in R. Use explicit local roots and safe argument passing. Check
readability before computation. Write only under the preparation output root.
Add fixture tests for substantive filtering and counting rules. Compare locally
recomputed results with current saved outputs wherever inputs are available.
Keep the original local checkout and paper files unchanged. Open a follow-up PR after the merge of PR 23 and check
GitHub Actions. Do not merge main or submit cloud jobs.
