# %% Current Figure 5: inputs and scope
#| Read the ancestry counts, ClinVar summaries, selected variants, and molecular-effect tables for the five panels.
stopifnot(requireNamespace("cowplot", quietly = TRUE))
# Use a font-aware temporary device when cowplot measures the Arial labels.
fig5_previous_null_device <- cowplot::set_null_device("png")
fig5_base <- "output/figure5_combined_draft"
fig5_print <- paste0(fig5_base, "/science_print_size")
fig5_rem <- paste0(fig5_print, "/stacked_left_arrows_REM")
fig5_populations <- c("afr", "amr", "eas", "eur", "mid", "sas")
fig5_inheritance <- c(AD = "Autosomal dominant", AR = "Autosomal recessive")
fig5_population_palette <- c(
  afr = "#CFDF99", amr = "#8ACE90", eas = "#F7B181",
  eur = "#47BB95", mid = "#F48C81", sas = "#F2D594"
)
fig5_class_palette <- c(
  "P/LP" = "#AF7AA1", VUS = "#F28E2C", Unclassified = "#4E79A7",
  "B/LB" = "#76B7B2"
)
fig5_classes <- names(fig5_class_palette)
fig5_class_shapes <- c("P/LP" = 16, VUS = 17, Unclassified = 15, "B/LB" = 18)

fig5_small_theme <- theme_classic(base_size = 6.5, base_family = "Arial") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 7, margin = margin(r = 3)),
    axis.text = element_text(size = 6.2, colour = "#333333"),
    axis.text.x = element_text(margin = margin(t = 2)),
    axis.line = element_line(linewidth = 0.22),
    axis.line.x = element_blank(),
    axis.ticks = element_line(linewidth = 0.18),
    axis.ticks.x = element_blank(),
    axis.ticks.length = grid::unit(0.6, "mm"),
    strip.text = element_text(size = 6.5, margin = margin(b = 3)),
    strip.background = element_blank(),
    plot.margin = margin(1, 1, 1, 1),
    legend.position = "none"
  )

# %% Panel A: rare variants and ClinVar classes
#| These tables contain the accepted counts across expression, splicing, and
#| protein QTLs. Each variant is assigned to the population with its highest
#| allele frequency. Counts are distinct within each inheritance group.
#| A gene or variant can occur in both inheritance groups. No new PIP filter
#| is applied here. The class bars use counts pooled across populations.
#| Class precedence in build_panel_A_v21.R is P/LP, then VUS, then B/LB.
#| A variant with no ClinVar classification is unannotated.
fig5_counts <- read_tsv(
  paper_file(paste0(fig5_print, "/data/A_counts.tsv")),
  show_col_types = FALSE
) %>%
  mutate(population = factor(population, levels = fig5_populations))
fig5_proportions <- read_tsv(
  paper_file(paste0(fig5_print, "/data/A_clinvar_proportions.tsv")),
  show_col_types = FALSE
) %>%
  mutate(clinvar_class = factor(
    clinvar_class, levels = c("P/LP", "VUS", "B/LB", "Unannotated")
  ))
fig5_mode_totals <- fig5_counts %>%
  group_by(allelicRequirements) %>%
  summarise(total = sum(total), .groups = "drop")
stopifnot(
  nrow(fig5_counts) == 12L,
  identical(fig5_mode_totals$total, c(1480, 3074)),
  all(fig5_proportions$total == if_else(
    fig5_proportions$allelicRequirements == "AD", 1480, 3074
  )),
  all(abs(fig5_proportions$proportion -
    fig5_proportions$n / fig5_proportions$total) < 1e-12)
)

fig5_A_counts <- ggplot(fig5_counts, aes(population, total, fill = population)) +
  geom_col(width = 0.85) +
  facet_grid(~allelicRequirements,
    labeller = labeller(allelicRequirements = fig5_inheritance)) +
  scale_fill_manual(values = fig5_population_palette, guide = "none") +
  scale_x_discrete(drop = FALSE, labels = str_to_upper) +
  scale_y_continuous(
    limits = c(0, 2000), breaks = seq(0, 2000, 500),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(y = "Fine-mapped rare variants") + fig5_small_theme +
  theme(
    panel.spacing.x = grid::unit(3, "mm"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
      size = 6.2, margin = margin(t = 2))
  )
fig5_A_classes <- ggplot(fig5_proportions) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1,
    fill = clinvar_class), colour = NA) +
  geom_text(
    data = fig5_proportions %>% filter(clinvar_class == "Unannotated"),
    aes(x = (xmin + xmax) / 2, y = 0.5, label = label),
    family = "Arial", size = 6.2 / .pt, colour = "white"
  ) +
  facet_grid(~allelicRequirements) +
  scale_fill_manual(values = c(
    fig5_class_palette[c("P/LP", "VUS", "B/LB")], Unannotated = "#4E79A7"
  ), drop = FALSE) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
    labels = scales::label_percent(accuracy = 1), expand = expansion(mult = 0)) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0)) +
  labs(x = "ClinVar classification (%)", y = NULL, fill = NULL) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_classic(base_size = 6.2, base_family = "Arial") +
  theme(
    axis.text.y = element_blank(), axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.text.x = element_text(size = 6.2, margin = margin(t = 2)),
    axis.title.x = element_text(size = 6.5, margin = margin(t = 3)),
    strip.text = element_blank(), strip.background = element_blank(),
    panel.spacing.x = grid::unit(3, "mm"), legend.position = "bottom",
    legend.text = element_text(size = 6.2),
    legend.key.width = grid::unit(2.1, "mm"),
    legend.key.height = grid::unit(2.1, "mm"),
    legend.key.spacing.x = grid::unit(0.8, "mm"),
    legend.box.spacing = grid::unit(1, "pt"),
    legend.margin = margin(1, 0, 0, 0), plot.margin = margin(0, 1, 0, 1)
  )
fig5_A <- fig5_A_counts / fig5_A_classes +
  plot_layout(heights = grid::unit(c(1, 3.5), c("null", "mm")))
# Keep the 0% and 100% labels inside their own facets, as in the archived SVG.
fig5_align_percent_labels <- function(grob) {
  if (inherits(grob, "text") &&
      identical(as.character(grob$label), c("0%", "50%", "100%"))) {
    grob$hjust <- c(0, 0.5, 1)
  }
  if (inherits(grob, "gtable")) {
    grob$grobs <- lapply(grob$grobs, fig5_align_percent_labels)
  }
  if (inherits(grob, "gTree") && length(grob$children)) {
    for (i in seq_along(grob$children)) {
      grob$children[[i]] <- fig5_align_percent_labels(grob$children[[i]])
    }
  }
  grob
}
fig5_A_grob <- fig5_align_percent_labels(patchwork::patchworkGrob(fig5_A))
fig5_A <- cowplot::ggdraw(fig5_A_grob)
save_plot(fig5_A, "Fig5_A_rebuilt", 72, 84)
print(fig5_A)
message("Panel A: AD = 1,480 and AR = 3,074 distinct rare variants.")

# %% Panel B: ClinVar coverage by population
#| This panel uses a different population assignment from panel A. A variant
#| is assigned to EUR if its EUR allele frequency is at least 75% of the
#| maximum frequency across the six populations. Other variants use the
#| population with the highest frequency. The overall AD and AR totals match
#| panel A, but some population totals differ. The dashed lines mark EUR.
#| The code also checks the saved two-sided Fisher tests for pooled non-EUR
#| versus EUR. The figure itself has no test labels or confidence intervals.
fig5_coverage <- read_tsv(
  paper_file(paste0(fig5_print, "/data/B_coverage.tsv")),
  show_col_types = FALSE
) %>%
  mutate(population = factor(population, levels = fig5_populations),
    label = sprintf("%.1f", percent))
fig5_B_tests_saved <- read_tsv(
  paper_file(paste0(fig5_base,
    "/panel_B_fisher_tests/pooled_non_EUR_vs_EUR.tsv")),
  show_col_types = FALSE
)
fig5_B_tests <- fig5_coverage %>%
  mutate(group = if_else(population == "eur", "EUR", "non_EUR")) %>%
  group_by(allelicRequirements, group) %>%
  summarise(annotated = sum(in_clinvar), unannotated = sum(total - in_clinvar),
    .groups = "drop") %>%
  group_by(allelicRequirements) %>%
  summarise(p_two_sided_rebuilt = fisher.test(cbind(
    annotated, unannotated
  ))$p.value, .groups = "drop") %>%
  left_join(fig5_B_tests_saved %>% select(allelicRequirements, p_two_sided),
    by = "allelicRequirements")
stopifnot(
  nrow(fig5_coverage) == 12L,
  all(abs(fig5_coverage$percent -
    100 * fig5_coverage$in_clinvar / fig5_coverage$total) < 1e-10),
  identical(fig5_coverage %>% group_by(allelicRequirements) %>%
    summarise(total = sum(total), .groups = "drop") %>% pull(total), c(1480, 3074)),
  all(abs(log(fig5_B_tests$p_two_sided_rebuilt /
    fig5_B_tests$p_two_sided)) < 1e-8)
)
fig5_B <- ggplot(fig5_coverage, aes(population, percent)) +
  geom_hline(data = fig5_coverage %>% filter(population == "eur"),
    aes(yintercept = percent), colour = fig5_population_palette[["eur"]],
    linetype = "dashed", linewidth = 0.25) +
  geom_segment(aes(xend = population, y = 0, yend = percent),
    colour = "#CDD1D3", linewidth = 0.3) +
  geom_point(aes(colour = emphasis), size = 1.5) +
  geom_text(aes(label = label, colour = emphasis), vjust = -1,
    size = 6.2 / .pt, family = "Arial") +
  facet_wrap(vars(allelicRequirements), ncol = 1, axes = "all_x",
    labeller = labeller(allelicRequirements = fig5_inheritance)) +
  scale_colour_manual(values = c(
    EUR = fig5_population_palette[["eur"]], `Other groups` = "#555555"
  ), guide = "none") +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 5),
    expand = expansion(mult = 0)) +
  scale_x_discrete(drop = FALSE, labels = str_to_upper) +
  labs(y = "In ClinVar (%)") + fig5_small_theme +
  theme(panel.spacing.y = grid::unit(6, "mm"))
save_plot(fig5_B, "Fig5_B_rebuilt", 45, 84)
print(fig5_B)
message("Panel B: both saved pooled Fisher tests match the aggregate counts.")

# %% Panel C: variant selection and pathogenic enrichment
#| The first three stages use the verified 2,928-variant starting pool.
#| Enrichment is the P/LP proportion at each stage divided by its proportion
#| in that pool. The last box counts nine AR genes. It requires a P/LP
#| reference and a distinct unclassified or VUS variant with high PIP.
#| This last search permits either effect direction. It is not a fourth
#| decreasing-effect filter. No new variant selection is run in this cell.
fig5_stages <- read_tsv(paper_file(
  "output/clingen_inheritance/ar_pathogenic_unclassified_funnel/stage_counts.tsv"
), show_col_types = FALSE)
fig5_verified_stages <- read_tsv(paper_file(paste0(
  fig5_base, "/stage_enrichment_tests/verified_stage_counts.tsv"
)), show_col_types = FALSE)
fig5_stage_tests_saved <- read_tsv(paper_file(paste0(
  fig5_base, "/stage_enrichment_tests/primary_enrichment_tests.tsv"
)), show_col_types = FALSE)
stopifnot(
  identical(fig5_stages$variants[1:3], c(2928, 526, 230)),
  identical(fig5_stages$genes, c(560, 375, 186, 9)),
  identical(fig5_stages$stage[1:3], fig5_verified_stages$stage),
  identical(fig5_stages$variants[1:3], fig5_verified_stages$variants),
  identical(fig5_stages$pathogenic_variants[1:3],
    fig5_verified_stages$pathogenic_variants)
)
fig5_stage_geometry <- tibble(
  stage = fig5_stages$stage, index = 1:4,
  ytop = c(10, 28.5, 47, 65.5), height = c(15, 15, 15, 17),
  top_width = c(43, 35, 28, 21), bottom_width = c(35, 28, 21, 14),
  colour = c("#DCECEF", "#76B5BE", "#155F6B", "#0A3943"),
  text_colour = c("#123E47", "#123E47", "white", "white")
) %>%
  left_join(fig5_stages, by = "stage") %>%
  mutate(enrichment_rebuilt = (pathogenic_variants / variants) /
    (first(pathogenic_variants) / first(variants)))
stopifnot(all(abs(fig5_stage_geometry$enrichment_rebuilt[1:3] -
  fig5_stages$enrichment[1:3]) < 1e-12))
fig5_stage_tests <- fig5_stage_geometry %>%
  filter(index %in% 2:3) %>%
  rowwise() %>%
  mutate(p_rebuilt = fisher.test(matrix(c(
    pathogenic_variants, variants - pathogenic_variants,
    21 - pathogenic_variants, 2928 - variants - (21 - pathogenic_variants)
  ), nrow = 2, byrow = TRUE))$p.value) %>%
  ungroup() %>%
  left_join(fig5_stage_tests_saved %>% select(selected_stage, p_value_two_sided),
    by = c("stage" = "selected_stage"))
stopifnot(all(abs(log(fig5_stage_tests$p_rebuilt /
  fig5_stage_tests$p_value_two_sided)) < 1e-8))
fig5_funnel_polygons <- fig5_stage_geometry %>%
  mutate(vertices = pmap(list(top_width, bottom_width, ytop, height),
    function(top_width, bottom_width, ytop, height) {
      tibble(x = 22 + c(-top_width, top_width, bottom_width, -bottom_width) / 2,
        y = ytop + c(0, 0, height, height))
    })) %>%
  select(index, colour, vertices) %>% unnest(vertices)
fig5_funnel_text <- fig5_stage_geometry %>% filter(index < 4) %>%
  mutate(filter_label = c("Rare eQTL + sQTL + pQTL", "PIP > 0.9",
    "Posterior mean < \u22120.3"),
    count_label = paste0(scales::comma(variants), " variants"),
    gene_label = paste0(genes, " AR genes"),
    enrichment_label = sprintf("%.1f\u00d7", enrichment_rebuilt))
fig5_C <- ggplot() +
  geom_polygon(data = fig5_funnel_polygons,
    aes(x, y, group = index, fill = colour), colour = NA) +
  geom_segment(data = fig5_funnel_text,
    aes(x = 22 + (top_width + bottom_width) / 4, xend = 51,
      y = ytop + 7.5, yend = ytop + 7.5),
    colour = "#DC6500", linewidth = 0.23) +
  geom_point(data = fig5_funnel_text, aes(x = 51, y = ytop + 7.5),
    colour = "#DC6500", size = 9.2) +
  geom_text(data = fig5_funnel_text,
    aes(x = 22, y = ytop + 3.5, label = filter_label, colour = text_colour),
    size = 6.2 / .pt, family = "Arial") +
  geom_text(data = fig5_funnel_text,
    aes(x = 22, y = ytop + 8.5, label = count_label, colour = text_colour),
    size = 8 / .pt, family = "Arial") +
  geom_text(data = fig5_funnel_text,
    aes(x = 22, y = ytop + 12.4, label = gene_label, colour = text_colour),
    size = 6.2 / .pt, family = "Arial") +
  geom_text(data = fig5_funnel_text,
    aes(x = 51, y = ytop + 7.5, label = enrichment_label),
    size = 7.5 / .pt, family = "Arial", colour = "white") +
  annotate("text", x = 51, y = 3.1,
    label = "Pathogenic\nvariant\nenrichment", size = 6.2 / .pt,
    family = "Arial", colour = "#123E47", lineheight = 0.95) +
  annotate("text", x = 22, y = 70.9,
    label = "Pathogenic +\nunclassified/VUS\nin the same gene",
    size = 6.2 / .pt, family = "Arial", colour = "white", lineheight = 0.95) +
  annotate("text", x = 22, y = 79.1,
    label = paste0(fig5_stages$genes[[4]], " AR genes"),
    size = 8 / .pt, family = "Arial", colour = "white") +
  annotate("rect", xmin = 10.5, xmax = 33.5, ymin = 64.5, ymax = 83.5,
    fill = NA, colour = "black", linewidth = 0.28) +
  scale_fill_identity() + scale_colour_identity() +
  scale_x_continuous(limits = c(0, 58), expand = expansion(mult = 0)) +
  scale_y_reverse(limits = c(84, 0), expand = expansion(mult = 0)) +
  coord_fixed(ratio = 1, clip = "off") + theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))
save_plot(fig5_C, "Fig5_C_rebuilt", 58, 84)
print(fig5_C)
message("Panel C: 2,928 -> 526 -> 230 variants; final endpoint = 9 AR genes.")

# %% Panels D and E: current REM effect displays
#| These tables contain the posterior mean QTL effects in the accepted plot.
#| Panel D has 14 variants at HEXA, PMM2, and GALM. Candidate effects have the
#| opposite sign to the pathogenic references. Panel E has 16 variants at six
#| other genes. Its candidate effects have the same sign as the references.
#| The plotted values are posterior effects, not log2 allelic fold changes.
#| Some benign comparators are common. PMM2 includes different assays, so the
#| effect magnitudes cannot be compared directly across all points.
#| The tables supply the current REM and ASJ [REM] labels. Effect direction
#| alone does not establish pathogenicity. No effect model is fitted here.
fig5_D_data <- read_tsv(paper_file(paste0(fig5_rem, "/data/D_display.tsv")),
  show_col_types = FALSE)
fig5_E_data <- read_tsv(paper_file(paste0(fig5_rem, "/data/E_display.tsv")),
  show_col_types = FALSE)
stopifnot(
  nrow(fig5_D_data) == 14L, nrow(fig5_E_data) == 16L,
  n_distinct(fig5_D_data$variant) == 14L,
  n_distinct(fig5_E_data$variant) == 16L,
  all(fig5_D_data$pip > 0.9), all(fig5_E_data$pip > 0.9),
  !any(str_detect(c(fig5_D_data$ancestry, fig5_E_data$ancestry), "OTH|ADMIX")),
  sum(str_detect(c(fig5_D_data$ancestry, fig5_E_data$ancestry), "REM")) == 4L
)
fig5_make_effect_panel <- function(data, gene_order, limits, breaks,
    gene_y, disease_y) {
  display <- data %>%
    mutate(symbol = factor(symbol, levels = gene_order),
      clinvar = factor(clinvar, levels = fig5_classes))
  groups <- display %>% group_by(symbol) %>%
    summarise(xmin = min(x) - 0.32, xmax = max(x) + 0.32,
      x = mean(x), .groups = "drop") %>%
    mutate(disease = recode(as.character(symbol),
      HEXA = "Tay Sachs", PMM2 = "Congenital disorder\nof glycosylation",
      GALM = "Galactosemia IV", GAA = "Pompe\nDisease",
      SMPD1 = "Acid sphingomyelinase\ndeficiency", CTH = "Cystathioninuria",
      PLG = "Hypoplasminogenemia", CC2D1A = "NDD",
      ACY1 = "Aminoacylase 1\ndeficiency"))
  ggplot(display, aes(x, effect)) +
    geom_segment(aes(xend = x, yend = 0), colour = "#C8CDD0", linewidth = 0.25) +
    geom_segment(data = groups, aes(x = xmin, xend = xmax, y = 0, yend = 0),
      inherit.aes = FALSE, colour = "#A8B0B5", linewidth = 0.32) +
    geom_point(aes(colour = clinvar, shape = clinvar), size = 2.2) +
    geom_text(aes(y = label_y, label = ancestry,
      vjust = if_else(effect > 0, 0, 1)), family = "Arial",
      size = 6.5 / .pt, colour = "#444444", lineheight = 0.95) +
    geom_text(data = groups, aes(x, y = gene_y, label = symbol),
      inherit.aes = FALSE, family = "Arial", size = 8 / .pt) +
    geom_text(data = groups, aes(x, y = disease_y, label = disease),
      inherit.aes = FALSE, family = "Arial", size = 6.5 / .pt,
      colour = "#444444", lineheight = 0.95) +
    scale_colour_manual(values = fig5_class_palette, drop = FALSE, guide = "none") +
    scale_shape_manual(values = fig5_class_shapes, drop = FALSE, guide = "none") +
    scale_x_continuous(expand = expansion(add = 0.72)) +
    scale_y_continuous(limits = limits, breaks = breaks,
      expand = expansion(mult = 0)) +
    labs(x = NULL, y = "QTL effect size") +
    theme_classic(base_size = 7, base_family = "Arial") +
    theme(
      axis.line.x = element_blank(), axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.line.y = element_line(linewidth = 0.22),
      axis.ticks.y = element_line(linewidth = 0.18),
      axis.ticks.length = grid::unit(0.7, "mm"),
      axis.text.y = element_text(size = 6.5, colour = "black"),
      axis.title.y = element_text(size = 7.5, margin = margin(r = 4)),
      plot.margin = margin(1, 5, 1, 1), legend.position = "none"
    )
}
fig5_D <- fig5_make_effect_panel(fig5_D_data, c("HEXA", "PMM2", "GALM"),
  c(-2.25, 1.68), seq(-2, 1, 0.5), 1.48, 1.19)
fig5_E <- fig5_make_effect_panel(fig5_E_data,
  c("GAA", "SMPD1", "CTH", "PLG", "CC2D1A", "ACY1"),
  c(-2.25, 0.64), seq(-2, 0, 0.5), 0.51, 0.26)
save_plot(fig5_D, "Fig5_D_rebuilt", 182, 51)
save_plot(fig5_E, "Fig5_E_rebuilt", 182, 51)
print(fig5_D)
print(fig5_E)
message("Panels D and E: all 30 effects and the four REM labels are retained.")

# %% Assemble: five-panel R reconstruction
#| The panel positions follow the accepted 184 x 225 mm stacked layout.
#| Black links leave the final funnel box and point left into D and E.
#| This R composition rebuilds all data panels. Font metrics, text positions,
#| endpoint label alignment, and the funnel frame can differ from the archived
fig5_width_mm <- 184
fig5_height_mm <- 225
fig5_place <- function(plot, x_mm, y_mm, width_mm, height_mm) {
  cowplot::draw_plot(plot,
    x = x_mm / fig5_width_mm,
    y = 1 - (y_mm + height_mm) / fig5_height_mm,
    width = width_mm / fig5_width_mm, height = height_mm / fig5_height_mm)
}
fig5_spine <- tribble(
  ~x, ~y, ~xend, ~yend,
  157.5, 80, 182, 80,
  182, 80, 182, 191.5
)
fig5_links <- tibble(x = 182, xend = 177.5,
  y = c(129.5, 191.5), yend = c(129.5, 191.5))
fig5_legend <- tibble(
  clinvar = factor(fig5_classes, levels = fig5_classes),
  x = c(48, 67, 86, 123), y = 221.2
)
fig5_rebuilt <- cowplot::ggdraw() +
  fig5_place(fig5_A_grob, 1, 6, 72, 84) +
  fig5_place(fig5_B, 76, 6, 45, 84) +
  fig5_place(fig5_C, 124, 6, 58, 84) +
  fig5_place(fig5_D, 1, 104, 182, 51) +
  fig5_place(fig5_E, 1, 166, 182, 51) +
  geom_segment(data = fig5_spine,
    aes(x = x / fig5_width_mm, xend = xend / fig5_width_mm,
      y = 1 - y / fig5_height_mm, yend = 1 - yend / fig5_height_mm),
    inherit.aes = FALSE, linewidth = 0.28) +
  geom_segment(data = fig5_links,
    aes(x = x / fig5_width_mm, xend = xend / fig5_width_mm,
      y = 1 - y / fig5_height_mm, yend = 1 - yend / fig5_height_mm),
    inherit.aes = FALSE, linewidth = 0.28,
    arrow = grid::arrow(length = grid::unit(1.8, "mm"), type = "closed")) +
  cowplot::draw_label("A", x = 1 / 184, y = 1 - 2.7 / 225,
    hjust = 0, size = 10, fontfamily = "Arial", fontface = "bold") +
  cowplot::draw_label("B", x = 76 / 184, y = 1 - 2.7 / 225,
    hjust = 0, size = 10, fontfamily = "Arial", fontface = "bold") +
  cowplot::draw_label("C", x = 124 / 184, y = 1 - 2.7 / 225,
    hjust = 0, size = 10, fontfamily = "Arial", fontface = "bold") +
  cowplot::draw_label("D", x = 1 / 184, y = 1 - 99.7 / 225,
    hjust = 0, size = 10, fontfamily = "Arial", fontface = "bold") +
  cowplot::draw_label("E", x = 1 / 184, y = 1 - 161.7 / 225,
    hjust = 0, size = 10, fontfamily = "Arial", fontface = "bold") +
  cowplot::draw_label("Opposite direction to pathogenic references",
    x = 92 / 184, y = 1 - 100 / 225, size = 7.5,
    fontfamily = "Arial", fontface = "bold") +
  cowplot::draw_label("Same direction as pathogenic references",
    x = 92 / 184, y = 1 - 162 / 225, size = 7.5,
    fontfamily = "Arial", fontface = "bold") +
  geom_point(data = fig5_legend,
    aes(x = x / fig5_width_mm, y = 1 - y / fig5_height_mm,
      colour = clinvar, shape = clinvar), inherit.aes = FALSE, size = 2.2) +
  geom_text(data = fig5_legend,
    aes(x = (x + 2.5) / fig5_width_mm, y = 1 - y / fig5_height_mm,
      label = clinvar), inherit.aes = FALSE, size = 7 / .pt,
    family = "Arial", hjust = 0) +
  scale_colour_manual(values = fig5_class_palette, guide = "none") +
  scale_shape_manual(values = fig5_class_shapes, guide = "none")
save_plot(fig5_rebuilt, "Fig5_rebuilt_layout", fig5_width_mm, fig5_height_mm)
options(repr.plot.width = 184 / 25.4, repr.plot.height = 225 / 25.4)
print(fig5_rebuilt)
cowplot::set_null_device(fig5_previous_null_device)

# %% Reference figure
#| Complete figure layout.
show_reference(5)
