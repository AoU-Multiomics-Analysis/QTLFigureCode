# Panel E — reserved slot. The bottom row of Figure 3 leaves space for a panel
# that has not been decided yet (see fig3/NOTES.md). Rendering a blank keeps the
# layout stable and lets assemble.R stay unconditional.
source("fig3/00_setup.R")

pE <- blank_plot

saveRDS(pE, file.path(OBJ_DIR, "E.rds"))
message("Panel E built (reserved slot).")
