# =============================================================================
# Figure 1 — end-to-end runner. From project root:  Rscript fig1/run_all.R
#   1. panels  : each stages its own inputs (fig1/inputs.tsv) and builds
#                -> fig1/objects/<panel>.rds
#   2. assemble: compose + save Figure1.pdf / Figure1.png
# Every panel loads ONLY the inputs it needs. To rebuild a single panel, run just
# that panel script (it re-stages/reads only its inputs) then fig1/assemble.R.
# =============================================================================
for (p in LETTERS[1:8]) source(sprintf("fig1/panels/panel%s.R", p))
source("fig1/assemble.R")
message("run_all complete.")
