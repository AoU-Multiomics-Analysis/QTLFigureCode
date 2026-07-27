# Panel C — conditional-TWAS QQ, colocalized vs matched non-colocalized (blood)
# Inputs: conditional_twas
source("fig4/00_setup.R")

qq_blood <- calculate_qq_data(qq_match_source(), traits = BLOOD_TRAITS)

p_blood_qq <- qq_blood %>%
  mutate(group = case_when(str_detect(group, "Matched") ~ "Null", TRUE ~ group)) %>%
  plot_qq() +
  theme(legend.position = "none")

saveRDS(p_blood_qq, file.path(OBJ_DIR, "C.rds"))
message("Panel C built.")
