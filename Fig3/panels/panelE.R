# Panel E — among colocalizing eQTL-pQTL variants, which annotations mark rarer
# and which mark commoner variants.
#
# Uses the SYMMETRIC colocalization background (fig3/symmetric.R): FastENLOC
# emits credible sets from one layer only, so the pQTL-specific side is
# reconstructed from the SuSiE table before anything is compared.
# Inputs: enloc_cs, susie_eqtl, susie_sqtl, susie_pqtl, gene_gtf
source("fig3/00_setup.R")

pE <- build_panelE_maf(coloc_maf_enrichment())

saveRDS(pE, file.path(OBJ_DIR, "E.rds"))
message("Panel E built.")
