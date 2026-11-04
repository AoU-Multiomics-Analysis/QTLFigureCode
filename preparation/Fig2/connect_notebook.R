source("preparation/common.R")
source("preparation/Fig2/methods.R")
source("preparation/Fig2/current_methods.R")
source("preparation/Fig3/plot_bridge.R")

plots <- readRDS(prep_reference("paper", "output/figure2_binned_draft/science_print_size/Science_display_plots.rds"))
cache <- set_names(LETTERS[1:6]) %>% map(~prepared_rds("source", paste0("Fig2/data/", .x, ".rds")))
frequency <- prepared_tsv("paper", "output/frequency_bins/frequency_bin_counts.tsv")
plots$A$data <- replace_prepared_frame(plots$A$data, frequency, "Fig2A frequency bins")
plots$B$data <- replace_prepared_frame(plots$B$data, cache$B$skew_df, "Fig2B unique variants")
for (entry in list(c(5, "skew_counts"), c(6, "pair_line_df"), c(7, "label_df"))) {
  i <- as.integer(entry[1]); plots$B$layers[[i]]$data <- replace_prepared_frame(plots$B$layers[[i]]$data, cache$B[[entry[2]]], paste0("Fig2B layer ", i))
}
plots$C$data <- replace_prepared_frame(plots$C$data, cache$C$afc_plot_df, "Fig2C annotation/aFC rows")
regression <- cache$C$log2_afc_regression
regression$x <- plots$C$layers[[3]]$data$x  # accepted annotation position, not an estimate
plots$C$layers[[3]]$data <- replace_prepared_frame(plots$C$layers[[3]]$data, regression, "Fig2C refitted regression annotations")
display <- fig2_cumulative_display(cache$D)
plots$D$data <- replace_prepared_frame(plots$D$data, display$plot_df, "Fig2D cumulative gene steps")
labels <- display$label_df %>% mutate(label = scales::percent(prop, accuracy = 1), y = if_else(direction == "Expression\ndecrease", 1.1, 4.65))
plots$D$layers[[5]]$data <- replace_prepared_frame(plots$D$layers[[5]]$data, labels, "Fig2D recomputed proportions")
plots$E$data <- replace_prepared_frame(plots$E$data, filter(cache$E, focal_lower_bound != -10), "Fig2E recalculated OR, CI, Fisher P")
plots$F$data <- replace_prepared_frame(plots$F$data, cache$F, "Fig2F protein-coding quantiles")
for (i in c(5, 6, 8, 10)) plots$F$layers[[i]]$data <- replace_prepared_subset(
  plots$F$layers[[i]]$data, cache$F, c("pct_bin", "gene_type"), paste0("Fig2F layer ", i))
prep_save_rds(plots, "paper", "output/figure2_binned_draft/science_print_size/Science_display_plots.rds")

# S16 uses newly derived ancestry bins and rolling summaries. The reference
# keeps only the accepted visual encoding, including the OTH-to-REM label.
bins <- prepared_rds("source", "Fig2/data/S_ancestry_stacked.rds")
relative <- "Supplemental Figures/01_Supplemental_drafts/Figure2_ancestry_enrichment_REM/Supplemental_Figure_S16_ancestry_enrichment_REM.rds"
supplement <- readRDS(prep_reference("paper", relative))
stack <- bins$stack_df %>% mutate(pop = toupper(gvs_max_subpop), pop = if_else(pop == "ADMIX", "OTH", pop))
supplement$data <- replace_prepared_frame(supplement$data, stack, "S16 ancestry bins")
supplement$layers[[2]]$data <- replace_prepared_frame(supplement$layers[[2]]$data, bins$line_df, "S16 rolling shared-frequency curve")
prep_save_rds(supplement, "paper", relative)
