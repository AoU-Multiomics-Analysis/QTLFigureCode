# Panel B — ancestry counts in cohort
# Inputs: comb_sample_list, ancestry_preds
source("fig1/00_setup.R")

sample_ids <- fread(input_path("comb_sample_list"))$research_id
ancestry_counts <- LoadAncestryCounts(input_path("ancestry_preds"), sample_ids)

pB <- PlotAncestryCounts(ancestry_counts)

saveRDS(pB, file.path(OBJ_DIR, "B.rds"))
message("Panel B built.")
