# Panel E — feature enrichment of colocalizing vs layer-specific variants,
# contrasted between assay pairs. See docs/superpowers/specs/
# 2026-07-27-fig3-panelE-coloc-enrichment-design.md
# Inputs: enloc_variant, susie_eqtl, susie_sqtl, susie_pqtl
source("fig3/00_setup.R")

pE <- build_panelE(coloc_feature_enrichment())

saveRDS(pE, file.path(OBJ_DIR, "E.rds"))
message("Panel E built.")
