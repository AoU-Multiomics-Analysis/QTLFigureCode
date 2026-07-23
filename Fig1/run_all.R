# =============================================================================
# Figure 1 — end-to-end runner. From project root:  Rscript fig1/run_all.R
#   1. prep    : read raw inputs once -> fig1/data/<panel>.rds
#   2. panels  : build each panel     -> fig1/objects/<panel>.rds
#   3. assemble: compose + save Figure1.pdf / Figure1.png
# To rebuild a single panel after editing it, just re-run that panel script then
# fig1/assemble.R — no need to re-read the raw data.
# =============================================================================
source("fig1/01_prep_data.R")
for (p in LETTERS[1:8]) source(sprintf("fig1/panels/panel%s.R", p))
source("fig1/assemble.R")
message("run_all complete.")
