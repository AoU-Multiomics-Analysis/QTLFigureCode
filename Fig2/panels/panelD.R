# Panel D — expected heterozygous expression from cumulative aFC (expression)
source("fig2/00_setup.R")
cumulative_effects <- readRDS(file.path(DATA_DIR, "D.rds"))

pD <- PlotExpectedExpression(cumulative_effects)

saveRDS(pD, file.path(OBJ_DIR, "D.rds"))
message("Panel D built.")
