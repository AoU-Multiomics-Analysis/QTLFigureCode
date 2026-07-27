# Fig5 end-to-end runner. Build every panel, then assemble. Run from project root:
#   Rscript fig5/run_all.R
for (p in LETTERS[1:7]) source(sprintf("fig5/panels/panel%s.R", p))
source("fig5/assemble.R")
message("run_all complete.")
