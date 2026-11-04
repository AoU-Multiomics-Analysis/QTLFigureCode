# Current QTL figure notebooks

The user approved the November 2, 2026 submission figures as the reference.
Update the five main R notebooks and add supplement notebooks when the saved
local inputs are complete. Preserve all existing edits in the original checkout.

Extend the existing notebook generator. Put current figure scripts and an input
index under `notebooks/current/`. Each notebook must contain executable R code,
explain its input files, and show the current figure. The code must distinguish
rebuilding a plot from repeating its upstream analysis. Do not claim a saved
plot is a new analysis. Preserve the original figure pipelines as historical code.

The local paper directory is supplied with `QTL_PAPER_ROOT`. The original figure
checkout is supplied with `QTL_SOURCE_ROOT`. Do not embed personal absolute paths.
Read controlled inputs in place. Write new outputs only under an ignored output
directory. Do not commit input tables or participant records. Commit aggregate
figure images in executed notebooks, as the existing repository permits.

Main references are the exact figures in the November 2 submission package.
Use the current post-QC seven-category Figure 1, compact Figure 2, current Figure
3, the Science Figure 4 with readable trait labels, and the Figure 5 REM version.
Add S4, S11-S16, S20-S25, S28, S30 and S31 from the source-data audit. Record the
missing and partial sections. S27 must retain a review note about its controls.

An alternative is to rewrite the six original data pipelines. That would change
more code than this notebook update needs. Another alternative is to replace
the notebooks with images alone. That would lose the executable figure code.
The selected approach preserves the current code and adds repeatable plot steps.

Check notebook generation, R syntax, source paths, aggregate-only outputs,
current section numbering, rendered outputs, and source-file hashes. Compare
the displayed main figures with the exact submission images. Use a review branch
and pull request. Do not merge into main as part of this task.
