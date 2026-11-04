source("preparation/common.R")
source("preparation/Fig3/methods.R")
source("preparation/Fig3/symmetric_methods.R")
source("preparation/Fig3/current_methods.R")
source("preparation/Fig3/plot_bridge.R")

plots <- set_names(LETTERS[1:5]) %>% map(~readRDS(prep_reference("source", paste0("Fig3/objects/", .x, ".rds"))))
panels <- set_names(LETTERS[1:4]) %>% map(~prepared_rds("source", paste0("Fig3/data/current_panel_", .x, ".rds")))
plots$A$data <- replace_prepared_frame(plots$A$data, panels$A, "Fig3A gene sharing")
shared <- sum(panels$A$prop[panels$A$NumColocs > 0])
shared_label <- sprintf("%.1f%% shared across ≥1 molecular layer", 100 * shared)
stopifnot(identical(plots$A$layers[[3]]$aes_params$label, shared_label))
plots$A$layers[[3]]$aes_params$label <- shared_label
full <- attr(plots$B, "grobs")$full
full$data <- replace_prepared_frame(full$data, panels$B$coloc_rate, "Fig3B per-assay credible-set rates")
full$patches$plots[[1]]$data <- replace_prepared_frame(full$patches$plots[[1]]$data, panels$B$signal_counts, "Fig3B component counts")
attr(plots$B, "grobs")$full <- full
plots$C$data <- replace_prepared_frame(plots$C$data, filter(panels$C, str_detect(term, ":")), "Fig3C refitted interactions")
plots$D$data <- replace_prepared_frame(plots$D$data, panels$D$data, "Fig3D transformed posterior effects")
slopes <- panels$D$slopes %>% mutate(label = paste0("β = ", round(slope, 2), "\nN = ", scales::comma(n)))
plots$D$layers[[6]]$data <- replace_prepared_frame(plots$D$layers[[6]]$data, slopes, "Fig3D refitted slopes")
enrichment <- prepared_rds("source", "Fig3/data/coloc_maf_enrichment.rds") %>% as_tibble() %>%
  filter(padj < 0.05, n >= PANEL_E_MIN_N, feature %in% names(PANEL_E_LABELS)) %>%
  mutate(dir = if_else(fold > 1, "rare", "common")) %>% arrange(fold) %>%
  mutate(label = factor(PANEL_E_LABELS[feature], levels = PANEL_E_LABELS[feature]))
plots$E$data <- replace_prepared_frame(plots$E$data, enrichment, "Fig3E refitted annotation effects")
iwalk(plots, ~prep_save_rds(.x, "source", paste0("Fig3/objects/", .y, ".rds")))
