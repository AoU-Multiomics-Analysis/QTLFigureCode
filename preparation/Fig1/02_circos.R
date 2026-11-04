# Prepare the exact reconciled circos data layers from the raw COMB tables.
# Geometry and label selection follow the archived current build_circos.R.
# The final SVG-only label-gap edits are a separate display composition step.
source("preparation/common.R")
config <- prep_init()
current <- "output/qtl_circos_simplified_options/N_main_figure_reconciled_20261027"
source_dir <- config$input_root
style <- "M"
width_mm <- 183
opening_degrees <- 45
chromosome_gap_degrees <- 0.7
fdr <- 0.05
retained_gene_types <- c("protein_coding", "lncRNA", "miRNA", "misc_RNA", "ribozyme", "scaRNA", "snRNA", "snoRNA")
pt_per_mm <- 72.27 / 25.4
# Match the named assay mapping in the existing Figure 1 panel H.
palette <- c(eQTL = "#00798C", sQTL = "#EDAE49", pQTL = "#D1495B")
background_size <- 0.24
background_alpha <- 0.42
simple_cap <- TRUE
cap_background_shape <- if (simple_cap) 16 else 17
cap_background_size <- if (simple_cap) background_size else 0.4
tracks <- tribble(
  ~assay, ~radius, ~height, ~cap, ~source,
  "eQTL", 52, 13, 200, "Fig1Data/eQTL/COMB.cis_qtl.txt.gz",
  "sQTL", 37, 13, 200, "Fig1Data/sQTL/COMB.cis_qtl.txt.gz",
  "pQTL", 22, 13, 150, "Fig1Data/pQTL/AoU_10750_COMB_pQTL.cis_qtl.txt.gz"
)


# Validate every raw file before reading it. Reuse source q values and coordinates.
walk(tracks$source, ~prep_input("input", .x))

# GRCh38 primary chromosome lengths. The local GENCODE v48 header specifies GRCh38.
chromosomes <- tibble(
  chromosome = paste0("chr", 1:22),
  length = c(248956422, 242193529, 198295559, 190214555, 181538259,
             170805979, 159345973, 145138636, 138394717, 133797422,
             135086622, 133275309, 114364328, 107043718, 101991189,
             90338345, 83257441, 80373285, 58617616, 64444167,
             46709983, 50818468)
) %>%
  mutate(
    chromosome_index = row_number(),
    span = length / sum(length) * (360 - opening_degrees - 21 * chromosome_gap_degrees),
    theta_start = 90 - opening_degrees / 2 -
      lag(cumsum(span + chromosome_gap_degrees), default = 0),
    theta_end = theta_start - span,
    theta_mid = (theta_start + theta_end) / 2
  )
gene_reference_path <- prep_input("paper", "output/qtl_circos_simplified_options/gene_reference.tsv")
gene_reference <- read_tsv(gene_reference_path, show_col_types = FALSE) %>%
  filter(chromosome %in% chromosomes$chromosome) %>%
  rename(gene_chromosome = chromosome, gene_start = start, gene_end = end)
stopifnot(!anyDuplicated(gene_reference$gene_id))

message("Read the combined-cohort lead cis-QTL summaries.")
points <- map2_dfr(tracks$assay, tracks$source, function(assay, source) {
  path <- file.path(source_dir, source)
  if (!file.exists(path)) stop("Missing input: ", path)
  # Preserve source text. The installed readr parser changes exponents below
  # approximately -302; base as.double() correctly reads these small values.
  read_tsv(path, show_col_types = FALSE,
           col_types = cols(pval_nominal = col_character(),
                            pval_beta = col_character(), qval = col_character())) %>%
    transmute(
      assay, phenotype_id, variant_id,
      gene_id = str_extract(phenotype_id, "ENSG[0-9]+"),
      chromosome = str_extract(variant_id, "^chr[0-9XY]+"),
      position = as.double(str_match(variant_id, "^chr[^:]+:([0-9]+)")[, 2]),
      source_pval_nominal = pval_nominal,
      source_pval_beta = pval_beta,
      source_qval = qval,
      pval_nominal = as.double(pval_nominal),
      pval_beta = as.double(pval_beta), qval = as.double(qval)
    )
}) %>%
  left_join(gene_reference, by = "gene_id", relationship = "many-to-one") %>%
  left_join(chromosomes, by = "chromosome", relationship = "many-to-one") %>%
  left_join(tracks %>% select(assay, radius, height, cap), by = "assay") %>%
  mutate(
    significant = qval < fdr,
    nominal_p_is_zero = pval_nominal == 0,
    neg_log10_p = -log10(pval_nominal),
    capped = neg_log10_p >= cap,
    plot_score = pmin(neg_log10_p, cap),
    theta = theta_start - position / length * span,
    plot_radius = radius + plot_score / cap * height,
    x = plot_radius * cos(theta * pi / 180),
    y = plot_radius * sin(theta * pi / 180),
    assay = factor(assay, levels = tracks$assay)
  )
stopifnot(
  nrow(points) == sum(map_int(tracks$source, ~ nrow(read_tsv(
    file.path(source_dir, .x), col_select = phenotype_id, show_col_types = FALSE)))),
  !anyDuplicated(points[c("assay", "phenotype_id")]),
  !anyNA(points$gene_id), !anyNA(points$theta),
  !anyNA(points$pval_nominal), !anyNA(points$qval),
  all(points$pval_nominal >= 0 & points$pval_nominal <= 1),
  all(points$qval >= 0 & points$qval <= 1),
  all(points$position > 0 & points$position <= points$length),
  all(points$chromosome == points$gene_chromosome, na.rm = TRUE),
  all(is.finite(points$x)), all(is.finite(points$y))
)

message("Retain GENCODE coding and RNA biotypes; exclude artifact annotations.")
biotype_audit <- points %>%
  mutate(retained = gene_type %in% retained_gene_types,
         gene_type = replace_na(gene_type, "unmapped")) %>%
  count(assay, gene_type, retained, name = "source_traits")
filter_summary <- biotype_audit %>%
  group_by(assay) %>%
  summarise(
    original_traits = sum(source_traits),
    retained_traits = sum(source_traits[retained]),
    excluded_traits = sum(source_traits[!retained]),
    percent_removed = 100 * excluded_traits / original_traits,
    .groups = "drop"
  )
points <- points %>% filter(gene_type %in% retained_gene_types)
stopifnot(!anyNA(points$gene_type),
          all(points$gene_type %in% retained_gene_types),
          nrow(points) == sum(filter_summary$retained_traits))
write_tsv(biotype_audit, prep_output("paper", file.path(current, "biotype_filter_audit.tsv")))
write_tsv(filter_summary, prep_output("paper", file.path(current, "filter_summary.tsv")))

if (style == "M") {
  message("Restrict to the paper's novel gene calls for the same QTL type.")
  novelty <- read_tsv(prep_input("paper", file.path(current, "paper_novelty_by_trait.tsv")),
                       col_types = cols(.default = col_character()))
  stopifnot(!anyDuplicated(novelty[c("assay", "phenotype_id")]))
  novel_keys <- novelty %>% filter(novel_for_figure == "TRUE") %>% select(assay, phenotype_id)
  novel_selection <- points %>%
    left_join(novelty %>% select(assay, phenotype_id, novel_in_paper),
              by = c("assay", "phenotype_id"), relationship = "one-to-one") %>%
    count(assay, novel_in_paper, name = "traits")
  stopifnot(!anyNA(novel_selection$novel_in_paper))
  write_tsv(novel_selection, prep_output("paper", file.path(current, "novel_selection_audit.tsv")))
  points <- points %>% semi_join(novel_keys, by = c("assay", "phenotype_id"))
  stopifnot(nrow(points) == nrow(novel_keys), all(points$qval < fdr),
            !any(points$assay == "pQTL" & points$gene_type == "lncRNA"))
}

summary <- points %>%
  group_by(assay) %>%
  summarise(
    source_traits = n(), source_genes = n_distinct(gene_id),
    significant_traits = sum(significant),
    significant_genes = n_distinct(gene_id[significant]),
    p_values_stored_as_zero = sum(nominal_p_is_zero),
    points_at_cap = sum(capped), missing_gene_name = sum(is.na(gene_name)),
    .groups = "drop"
  )

# Labels are a separate design input. Filter them after the novel targets exist.
label_path <- prep_input("paper", "output/qtl_circos_simplified_options/L_noncoding_highlights/gene_labels.tsv")
requested_labels <- read_tsv(label_path, show_col_types = FALSE) %>%
  semi_join(points %>% distinct(genes = gene_name), by = "genes")
stopifnot(nrow(requested_labels) == 26)
write_tsv(requested_labels, prep_output("paper", file.path(current, "gene_labels.tsv")))
write_tsv(novelty %>% filter(novel_for_figure == "TRUE", gene_name %in% requested_labels$genes),
          prep_output("paper", file.path(current, "selected_gene_evidence.tsv")))
label_genes <- requested_labels %>% separate_longer_delim(genes, delim = ";")
label_points <- points %>%
  filter(significant) %>%
  inner_join(label_genes, by = c("gene_name" = "genes"), relationship = "many-to-one") %>%
  mutate(point_fill = unname(palette[as.character(assay)]))
stopifnot(setequal(label_points$label, requested_labels$label))
anchors <- label_points %>%
  arrange(assay, desc(plot_score), qval, phenotype_id) %>%
  group_by(label) %>% slice_head(n = 1) %>% ungroup() %>%
  arrange(theta) %>%
  mutate(
    label_theta = theta,
    anchor_assay = as.character(assay),
    anchor_variant = variant_id
  )

# Move nearby labels apart without changing their order around the genome.
if (nrow(anchors) > 1) {
  for (i in seq.int(2, nrow(anchors))) {
    anchors$label_theta[[i]] <- max(anchors$label_theta[[i]], anchors$label_theta[[i - 1]] + 3.3)
  }
  anchors$label_theta <- anchors$label_theta - mean(anchors$label_theta - anchors$theta)
}
anchors <- anchors %>% select(-label_theta) %>%
  left_join(read_tsv(prep_input("paper", "output/qtl_circos_simplified_options/L_noncoding_highlights/label_anchors.tsv"),
                     show_col_types = FALSE) %>% select(label, label_theta),
            by = "label", relationship = "one-to-one")
stopifnot(!anyNA(anchors$label_theta))
anchors <- anchors %>%
  mutate(
    sx = 72.0 * cos(theta * pi / 180), sy = 72.0 * sin(theta * pi / 180),
    ex = 78.2 * cos(label_theta * pi / 180), ey = 78.2 * sin(label_theta * pi / 180),
    tx = 79.1 * cos(label_theta * pi / 180), ty = 79.1 * sin(label_theta * pi / 180),
    text_angle = if_else(cos(label_theta * pi / 180) < 0, label_theta + 180, label_theta),
    text_hjust = if_else(cos(label_theta * pi / 180) < 0, 1, 0)
  )

# One dot per novel gene/QTL-type pair. Put the dots between the straight
# connector and the label, so the annotation does not extend the outer margin.
# The order from the circle toward the gene name is eQTL, sQTL, pQTL.
gene_type_dots <- label_points %>%
  group_by(label, assay) %>%
  summarise(gene_ids = paste(sort(unique(gene_id)), collapse = ";"), .groups = "drop") %>%
  arrange(label, assay) %>%
  group_by(label) %>%
  mutate(dot_index = row_number(), n_types = n()) %>% ungroup() %>%
  left_join(anchors %>% select(label, label_theta), by = "label", relationship = "many-to-one") %>%
  mutate(
    dot_radius = 77.25 - (n_types - dot_index) * 1.5,
    x = dot_radius * cos(label_theta * pi / 180),
    y = dot_radius * sin(label_theta * pi / 180),
    dot_fill = unname(palette[as.character(assay)])
  )
stopifnot(!anyDuplicated(gene_type_dots[c("label", "assay")]),
          nrow(gene_type_dots) == nrow(label_points %>% distinct(label, assay)),
          all(gene_type_dots$dot_radius > 71),
          all(abs(gene_type_dots$x) + 0.7 < 91.5),
          all(abs(gene_type_dots$y) + 0.7 < 91.5))
write_tsv(gene_type_dots, prep_output("paper", file.path(current, "gene_label_dot_types.tsv")))
write_tsv(summary, prep_output("paper", file.path(current, "summary.tsv")))
write_tsv(points, prep_output("paper", file.path(current, "plotted_lead_variants.tsv.gz")))
write_tsv(anchors %>% select(label, reason, gene_name, gene_id, chromosome,
                            anchor_assay, anchor_variant, position, theta, label_theta,
                            sx, sy, ex, ey, tx, ty, text_angle, text_hjust),
          prep_output("paper", file.path(current, "label_anchors.tsv")))
write_tsv(label_points %>% select(label, gene_name, gene_id, assay, phenotype_id,
                                 variant_id, chromosome, position, pval_nominal, qval, capped),
          prep_output("paper", file.path(current, "highlighted_gene_signals.tsv")))
axis_scales <- points %>% group_by(assay) %>%
  summarise(minimum = 0, maximum = first(cap), points = n(),
            capped_points = sum(capped), capped_percent = 100 * mean(capped),
            .groups = "drop")
write_tsv(axis_scales, prep_output("paper", file.path(current, "axis_scales.tsv")))

arc <- function(chr, inner, outer) {
  angles <- seq(chr$theta_start, chr$theta_end, length.out = 100) * pi / 180
  tibble(x = c(inner * cos(angles), outer * cos(rev(angles))),
         y = c(inner * sin(angles), outer * sin(rev(angles))))
}
bands <- map_dfr(seq_len(nrow(chromosomes)), function(i) {
  chr <- chromosomes %>% slice(i)
  arc(chr, 67.7, 71) %>% mutate(group = i, fill = if_else(i %% 2 == 0, "#E4E5E7", "#F0F1F2"))
})
backgrounds <- crossing(chromosomes, tracks %>% select(assay, radius, height)) %>%
  pmap_dfr(function(chromosome, length, chromosome_index, span, theta_start,
                   theta_end, theta_mid, assay, radius, height) {
    arc(tibble(theta_start, theta_end), radius, radius + height) %>%
      mutate(group = paste(assay, chromosome),
             fill = "#FFFFFF")
  })
gridlines <- crossing(chromosomes, tracks %>% select(assay, radius, height, cap)) %>%
  mutate(value = map(cap, ~ seq(0, .x, by = 50))) %>% unnest(value) %>%
  pmap_dfr(function(chromosome, length, chromosome_index, span, theta_start,
                   theta_end, theta_mid, assay, radius, height, cap, value) {
    theta <- seq(theta_start, theta_end, length.out = 100) * pi / 180
    r <- radius + value / cap * height
    tibble(x = r * cos(theta), y = r * sin(theta), group = paste(assay, chromosome, value))
  })
chromosome_text <- chromosomes %>%
  mutate(x = 69.35 * cos(theta_mid * pi / 180), y = 69.35 * sin(theta_mid * pi / 180),
         text = str_remove(chromosome, "chr"))
ticks <- chromosomes %>%
  mutate(position = map(length, ~ seq(0, .x, by = 1e8))) %>%
  unnest(position) %>%
  mutate(theta = theta_start - position / length * span,
         x = 71.1 * cos(theta * pi / 180), y = 71.1 * sin(theta * pi / 180),
         xend = 72.1 * cos(theta * pi / 180), yend = 72.1 * sin(theta * pi / 180),
         tx = 73.7 * cos(theta * pi / 180), ty = 73.7 * sin(theta * pi / 180),
         text_angle = theta - 90,
         text_angle = if_else(cos(theta * pi / 180) < 0, text_angle + 180, text_angle),
         text = as.character(position / 1e6))
# Follow the left edge of the opening. The small perpendicular offset keeps
# the axis and its tick labels in the white gap, clear of the data points.
axis_theta_degrees <- 90 + opening_degrees / 2
axis_theta <- axis_theta_degrees * pi / 180
axis_offset <- 0.8
normal_x <- sin(axis_theta)
normal_y <- -cos(axis_theta)
axis_positions <- tracks %>% select(assay, radius, height, cap) %>%
  mutate(theta_degrees = axis_theta_degrees, offset = axis_offset,
         sx = radius * cos(axis_theta) + axis_offset * normal_x,
         sy = radius * sin(axis_theta) + axis_offset * normal_y,
         ex = (radius + height) * cos(axis_theta) + axis_offset * normal_x,
         ey = (radius + height) * sin(axis_theta) + axis_offset * normal_y)
axes <- axis_positions %>%
  mutate(value = map(cap, ~ seq(0, .x, by = 50))) %>% unnest(value) %>%
  mutate(tick_radius = radius + height * value / cap,
         x = tick_radius * cos(axis_theta) + axis_offset * normal_x,
         y = tick_radius * sin(axis_theta) + axis_offset * normal_y,
         xend = x + normal_x, yend = y + normal_y,
         tx = x + 1.6 * normal_x, ty = y + 1.6 * normal_y,
         text_angle = axis_theta_degrees - 90)
write_tsv(axis_positions, prep_output("paper", file.path(current, "axis_positions.tsv")))
write_tsv(axes, prep_output("paper", file.path(current, "axis_ticks.tsv")))
leaders <- anchors %>% select(label, sx, sy, label_theta) %>%
  left_join(gene_type_dots %>% group_by(label) %>%
              summarise(end_radius = min(dot_radius) - 0.85, .groups = "drop"),
            by = "label", relationship = "one-to-one") %>%
  mutate(ex = end_radius * cos(label_theta * pi / 180),
         ey = end_radius * sin(label_theta * pi / 180)) %>%
  select(label, sx, sy, ex, ey)
write_tsv(leaders, prep_output("paper", file.path(current, "gene_connector_segments.tsv")))

legend_rows <- summary %>%
  mutate(y = 2.7 - (row_number() - 1) * 5.2,
         count_text = scales::comma(significant_traits))
legend_layout <- bind_rows(
  tibble(role = "header", label = "Novel \u00b7 q < 0.05", x = 0, y = 8.2),
  legend_rows %>% transmute(role = "assay", label = as.character(assay), x = -7.2, y),
  legend_rows %>% transmute(role = "count", label = count_text, x = 11, y),
  legend_rows %>% transmute(role = "swatch", label = as.character(assay), x = -10, y)
)
write_tsv(legend_layout, prep_output("paper", file.path(current, "legend_layout.tsv")))

# A panel-flow arrow begins just above the centered legend, passes through
# the opening, and turns right. There is no arrowhead at the legend end.
panel_connector <- tibble(vertex = 1:3, x = c(0, 0, 88), y = c(11.5, 86, 86))
write_tsv(panel_connector, prep_output("paper", file.path(current, "panel_connector_vertices.tsv")))


expected <- read_tsv(prep_reference("paper", file.path(current, "plotted_lead_variants.tsv.gz")),
                     col_types = cols(.default = col_character()), progress = FALSE)
keys <- c("assay", "phenotype_id", "variant_id")
check <- points %>% mutate(assay = as.character(assay)) %>%
  select(all_of(keys), plot_score, x, y) %>%
  full_join(expected %>% transmute(assay, phenotype_id, variant_id,
    plot_score = as.double(plot_score), x = as.double(x), y = as.double(y)),
    by = keys, suffix = c(".new", ".saved"), relationship = "one-to-one")
stopifnot(nrow(check) == nrow(points), !anyNA(check),
  max(abs(check$plot_score.new - check$plot_score.saved)) < 1e-12,
  max(abs(check$x.new - check$x.saved)) < 1e-12,
  max(abs(check$y.new - check$y.saved)) < 1e-12)
message("Figure 1 circos inputs verified: ", nrow(points), " lead points and 26 labels match the saved data.")
