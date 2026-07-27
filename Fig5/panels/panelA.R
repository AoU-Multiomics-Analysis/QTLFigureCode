# Fig5 Panel A — fine-mapped rare ClinGen variants by ancestry (AD/AR) over a
# per-inheritance ClinVar-classification proportion bar. Untagged; run from root.
source("fig5/00_setup.R")

ref  <- get_reference()
eqtl <- get_eqtl()

A <- build_panelA(eqtl, ref$clingen)
saveRDS(A, file.path(OBJ_DIR, "A.rds"))
message("Panel A built.")
