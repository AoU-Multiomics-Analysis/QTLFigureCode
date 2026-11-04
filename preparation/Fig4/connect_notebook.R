source("preparation/common.R")
source("preparation/Fig3/plot_bridge.R")

for (panel in c("B", "C", "D")) {
  plot <- readRDS(prep_reference("source", paste0("Fig4/objects/", panel, ".rds")))
  data <- prepared_rds("source", paste0("Fig4/data/current_panel_", panel, ".rds"))
  if (panel == "C") data <- data %>% mutate(group = if_else(str_detect(as.character(group), "Matched"), "Null", as.character(group)))
  if (panel == "D") data <- data$plot_data
  plot$data <- replace_prepared_frame(plot$data, data, paste0("Fig4", panel, " recalculated data/statistics"))
  prep_save_rds(plot, "source", paste0("Fig4/objects/", panel, ".rds"))
}

# A remains unavailable when its full harmonized gene result is missing.
# Do not copy the saved A object into the generated overlay in that case.
traits_path <- file.path(.prep_state$config$output_root, "source", "Fig4/data/current_panel_A.tsv")
if (file.exists(traits_path)) {
  traits <- prepared_tsv("source", "Fig4/data/current_panel_A.tsv")
  plot <- readRDS(prep_reference("source", "Fig4/objects/A_labeled.rds"))
  replacement <- plot$data %>% mutate(trait_raw = sub("___.*$", "", as.character(trait))) %>%
    left_join(traits %>% select(trait_raw = trait, trait_category, count, pct_coloc_union),
              by = c("trait_raw", "trait_category"), relationship = "many-to-one") %>%
    mutate(value = if_else(metric == "a_count", count, pct_coloc_union)) %>% select(all_of(names(plot$data)))
  plot$data <- replace_prepared_frame(plot$data, replacement, "Fig4A gene counts and lenient rates")
  for (i in seq_along(plot$layers)) {
    reference <- plot$layers[[i]]$data
    if (is.data.frame(reference) && all(c("trait_raw", "count") %in% names(reference))) {
      statistics <- c("count", "pct_stringent_union", "pct_lenient_union", "pct_coloc_union", "rate_type")
      replacement <- reference %>% select(-any_of(statistics)) %>%
        left_join(traits %>% select(trait_raw = trait, trait_category, any_of(statistics)),
                  by = c("trait_raw", "trait_category"), relationship = "many-to-one")
      plot$layers[[i]]$data <- replace_prepared_frame(reference, replacement, "Fig4A trait-label statistics")
    }
  }
  prep_save_rds(plot, "source", "Fig4/objects/A_labeled.rds")
} else message("Fig4A has no prepared overlay: its full harmonized gene input is missing")
