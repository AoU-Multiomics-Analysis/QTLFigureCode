# Panel D — novel eGene + sGene counts vs external consortia (pGene deferred)
source("fig1/00_setup.R")
d <- readRDS(file.path(DATA_DIR, "D.rds"))

pD <- PlotNovelQTLCount(d)

saveRDS(pD, file.path(OBJ_DIR, "D.rds"))
message("Panel D built.")
