# Fig5 Panel E — ACVRL1 causal promoter variants, PromoterAI base-resolution
# effect (signed-square y-axis). Untagged; run from the project root.
source("fig5/00_setup.R")

eqtl <- get_eqtl()
vf   <- get_variant_freq()

E <- build_panelE(eqtl, vf)
saveRDS(E, file.path(OBJ_DIR, "E.rds"))
message("Panel E built.")
