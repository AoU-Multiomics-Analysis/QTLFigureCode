# Fig5 Panel B — |average QTL effect| by variant-effect category (Control/AD/AR),
# with 95% CIs. Untagged; run from the project root.
source("fig5/00_setup.R")

ref  <- get_reference()
eqtl <- get_eqtl()

B <- build_panelB(eqtl, ref$clingen)
saveRDS(B, file.path(OBJ_DIR, "B.rds"))
message("Panel B built.")
