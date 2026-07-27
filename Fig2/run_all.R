# =============================================================================
# Figure 2 — end-to-end runner. From project root:  Rscript fig2/run_all.R
#   1. prep   : read raw inputs once -> fig2/data/<panel>.rds
#   2. panels : build each panel     -> fig2/objects/<panel>.rds
#   3. assemble: compose + save Figure2.pdf / Figure2.png
# To rebuild a single panel after editing it, just re-run that panel script then
# fig2/assemble.R — no need to re-read the raw data.
# =============================================================================
source("fig2/01_prep_data.R")
for (p in LETTERS[1:7]) source(sprintf("fig2/panels/panel%s.R", p))
source("fig2/assemble.R")
message("run_all complete.")
