# =============================================================================
# Supplemental Figure 1 — end-to-end runner. From the project root:
#   Rscript suppfig1/run_all.R
#   1. panels  : each stages the inputs it needs (suppfig1/inputs.tsv) and builds
#                the derived tables it needs (memoized in suppfig1/data/)
#   2. assemble: compose + save SupplementalFig1.pdf / .png
# The extra/ scripts are standalone: Rscript suppfig1/extra/<name>.R
# The derived-data cache is keyed by NAME: clear suppfig1/data/ after changing an input.
# =============================================================================
for (p in LETTERS[1:4]) source(sprintf("suppfig1/panels/panel%s.R", p))
source("suppfig1/assemble.R")
message("run_all complete.")
