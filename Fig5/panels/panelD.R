# Fig5 Panel D — ACVRL1 eQTL fine-mapping (regional PIP + zoom inset + gene track).
# Untagged; run from the project root.
source("fig5/00_setup.R")

eqtl <- get_eqtl()

D <- build_panelD(eqtl)
saveRDS(D, file.path(OBJ_DIR, "D.rds"))
message("Panel D built.")
