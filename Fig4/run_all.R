# =============================================================================
# Figure 4 — end-to-end runner. From the project root:  Rscript fig4/run_all.R
#   1. panels  : each stages the inputs it needs (fig4/inputs.tsv) and builds the
#                derived tables it needs (memoized in fig4/data/) -> fig4/objects/
#   2. assemble: compose + save Figure4.pdf / Figure4.png, and the labelled-A
#                variant Figure4_labeled.pdf / .png
# The derived-data cache is keyed by NAME, not content: clear fig4/data/ after
# changing an input file.
# =============================================================================
for (p in LETTERS[1:4]) source(sprintf("fig4/panels/panel%s.R", p))
source("fig4/panels/panelA_labeled.R")
source("fig4/assemble.R")
message("run_all complete.")
