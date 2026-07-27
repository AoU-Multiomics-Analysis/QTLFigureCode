# Fig5 Panel C — log2 aFC of the strongest decreasing variant per gene (top 7 per
# AD/AR). Uses the optional private no-homozygote list if present. Untagged.
source("fig5/00_setup.R")

ref  <- get_reference()
eqtl <- get_eqtl()

C <- build_panelC(eqtl, ref$clingen, nohom = load_nohomozygotes())
saveRDS(C, file.path(OBJ_DIR, "C.rds"))
message("Panel C built.")
