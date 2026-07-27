# Fig5 Panel F — SPR eQTL fine-mapping (regional PIP + zoom inset + gene track).
# Untagged; run from the project root.
source("fig5/00_setup.R")

eqtl <- get_eqtl()

F_panel <- build_panelF(eqtl)
saveRDS(F_panel, file.path(OBJ_DIR, "F.rds"))
message("Panel F built.")
