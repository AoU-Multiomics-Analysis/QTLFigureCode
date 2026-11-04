# %% Colocalization membership transformations
#| Read source memberships, remove duplicate gene–trait–assay rows, and exclude the internal union layer.
# @include preparation/Fig4/methods.R
# @include preparation/Fig4/current_methods.R

# %% Read source gene–trait memberships
#| QTL_FIG4_MEMBERSHIP can select the current membership result. The source table contains gene–trait records, not participant records.
membership_input <- prep_input_file("QTL_FIG4_MEMBERSHIP", "input", c(
  "high_level_coloc_outputs/final_summaries/colocalizing_genes_long.tsv",
  "Fig4Data/high_level_coloc_outputs/final_summaries/colocalizing_genes_long.tsv"))
coloc_membership <- read_any(membership_input)
upset <- fig4_upset_counts(coloc_membership)
intersections <- upset$intersections %>%
  arrange(desc(gene_trait_pairs)) %>%
  mutate(index = row_number(), sets = strsplit(qtl_combination, "+", fixed = TRUE),
         n_sets = lengths(sets), multilayer = n_sets >= 2)
stopifnot(nrow(intersections) == 7L, upset$multilayer_n == 1617L,
          sum(intersections$gene_trait_pairs[intersections$multilayer]) == 1617L)
qtl_sets <- c("eQTL", "pQTL", "sQTL")
message("Calculated seven exclusive QTL membership intersections.")

# %% Reconstruct the aggregate UpSet matrix
#| One column is one exclusive QTL combination. Filled dots identify its component assays. The top bars give the calculated exclusive counts. The side bars sum all combinations that contain an assay.
matrix <- tidyr::crossing(index = intersections$index, assay = qtl_sets) %>%
  left_join(intersections %>% select(index, sets), by = "index") %>%
  mutate(active = purrr::map2_lgl(assay, sets, ~.x %in% .y),
         y = match(assay, rev(qtl_sets)))
connectors <- matrix %>% filter(active) %>% group_by(index) %>%
  summarise(low = min(y), high = max(y), .groups = "drop")
set_sizes <- matrix %>% filter(active) %>%
  left_join(intersections %>% select(index, gene_trait_pairs), by = "index") %>%
  group_by(assay, y) %>% summarise(gene_trait_pairs = sum(gene_trait_pairs), .groups = "drop")
bars <- ggplot(intersections, aes(index, gene_trait_pairs, fill = multilayer)) +
  geom_col(width = 0.65) + geom_text(aes(label = scales::comma(gene_trait_pairs)), vjust = -0.35, size = 2.8) +
  scale_fill_manual(values = c(`FALSE` = "#2F5597", `TRUE` = "#CC79A7"), guide = "none") +
  scale_x_continuous(limits = c(0.4, 7.6), breaks = NULL, expand = expansion(mult = 0)) +
  scale_y_continuous(labels = scales::label_comma(), expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Gene–trait pairs per QTL combination") + theme_classic(base_size = 9) +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank())
dots <- ggplot(matrix, aes(index, y)) +
  geom_point(size = 2.5, colour = "#DDDDDD") +
  geom_segment(data = connectors, aes(x = index, xend = index, y = low, yend = high), inherit.aes = FALSE, linewidth = 0.7) +
  geom_point(data = filter(matrix, active), size = 2.5, colour = "#333333") +
  scale_x_continuous(limits = c(0.4, 7.6), breaks = NULL, expand = expansion(mult = 0)) +
  scale_y_continuous(limits = c(0.5, 3.5), breaks = 1:3, labels = rev(qtl_sets), expand = expansion(mult = 0)) +
  labs(x = NULL, y = NULL) + theme_void(base_size = 9) + theme(axis.text.y = element_text())
totals <- ggplot(set_sizes, aes(gene_trait_pairs, y, fill = assay)) +
  geom_col(orientation = "y", width = 0.55) +
  scale_fill_manual(values = c(eQTL = "#0072B2", pQTL = "#D55E00", sQTL = "#009E73"), guide = "none") +
  scale_x_reverse(labels = scales::label_comma()) +
  scale_y_continuous(limits = c(0.5, 3.5), breaks = NULL, expand = expansion(mult = 0)) +
  labs(x = "Gene–trait pairs", y = NULL) + theme_classic(base_size = 8) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())

# %% Display and save the aggregate intersections
#| Pink bars identify combinations with at least two QTL classes. The counts come from distinct source gene–trait memberships. Each pair contributes to one exclusive combination.
figure <- patchwork::plot_spacer() + bars + totals + dots +
  patchwork::plot_layout(ncol = 2, widths = c(1, 3.2), heights = c(3, 1.2))
print(figure)
save_plot(figure, "supplements/S24", 184, 145)
readr::write_tsv(intersections %>% select(qtl_combination, gene_trait_pairs, n_sets), output_file("supplements/S24_intersections.tsv"))
