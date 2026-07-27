# Fig5 Panel G — SPR Borzoi multi-track allelic browser + CTCF attribution logo.
# A cowplot object wrapped for patchwork. Untagged; run from the project root.
source("fig5/00_setup.R")

eqtl <- get_eqtl()
vf   <- get_variant_freq()

G <- build_panelG(eqtl, vf)
saveRDS(G, file.path(OBJ_DIR, "G.rds"))
message("Panel G built.")
