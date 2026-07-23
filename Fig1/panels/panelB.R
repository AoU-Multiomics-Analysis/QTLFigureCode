# Panel B — ancestry counts in cohort
source("fig1/00_setup.R")
d <- readRDS(file.path(DATA_DIR, "B.rds"))

pB <- PlotAncestryCounts(d)

saveRDS(pB, file.path(OBJ_DIR, "B.rds"))
message("Panel B built.")
