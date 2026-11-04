# Visible notebook data processing

The user wants the notebooks to follow the scripts and show data processing.
Each supported notebook must read the original local result tables, show its
filters, joins, counts and models, and construct its plots in R. A saved data
or plot RDS must not supply a notebook result. Source inputs remain outside Git.

Use the existing R transformations and current figure specifications. Expand
shared R files into visible notebook cells through explicit `# @include` lines
in the authored R scripts. Keep intermediate objects in the current R session.
The alternatives, linking preparation scripts or replacing data inside saved
plots, were already used and do not meet the user's clarified requirement.

Named local roots remain `QTL_INPUT_ROOT`, `QTL_SUSIE_ROOT`, `QTL_PAPER_ROOT`,
and `QTL_SOURCE_ROOT`. Validate and hash actual input files. Preserve scientific
selection rules and compare generated values with the accepted results locally.
Saved objects may be read by separate validation commands, but never by the
notebook or its plot builders. The reference PNG cells can remain.

Missing original inputs must be stated in the guide. Figure 4A can show its
processing code and skip that panel with a status message when the full input
is absent. Pre-QC availability and S30 retain their supplied aggregate-input
boundaries. Do not invent source calculations that cannot be recovered.

Use tidyverse, clean plots, no titles or subtitles, and no notebook Scope or
submission-version boilerplate. Keep outputs limited to aggregate images and
status messages. Do not change the original figure checkout or paper inputs.
Publish a new PR. Do not merge or submit scientific cloud jobs.
