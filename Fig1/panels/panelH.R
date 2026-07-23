# Panel H — fine-mapping functional enrichment by assay
source("fig1/00_setup.R")
d <- readRDS(file.path(DATA_DIR, "H.rds"))

pH <- PlotFinemappingEnrichment(d$enrichment, d$functional_groups)

saveRDS(pH, file.path(OBJ_DIR, "H.rds"))
message("Panel H built.")
