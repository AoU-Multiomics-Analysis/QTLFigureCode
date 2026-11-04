# %% Compute circos positions and selected gene labels
#| Positions come from the raw variant rows. GRCh38 chromosome lengths determine the sectors.
#| Gene selection and manual label angles are design inputs; point and connector coordinates are calculated here.
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
gene_reference <- genes %>% filter(chromosome %in% chromosomes$chromosome) %>%
  rename(gene_chromosome = chromosome, gene_start = start, gene_end = end)
# Derive angles and capped radial scores from the original nominal P value and variant position.
points <- raw %>% transmute(assay, phenotype_id, variant_id,
  gene_id = str_extract(phenotype_id, "ENSG[0-9]+"),
  chromosome = str_extract(variant_id, "^chr[0-9XY]+"),
  position = as.double(str_match(variant_id, "^chr[^:]+:([0-9]+)")[, 2]),
  source_pval_nominal, source_pval_beta, source_qval,
  pval_nominal = as.double(source_pval_nominal),
  pval_beta = as.double(source_pval_beta), qval = as.double(source_qval)) %>%
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
points <- points %>% filter(gene_type %in% retained_gene_types) %>%
  semi_join(novelty %>% filter(novel_for_figure), by = c("assay", "phenotype_id"))
stopifnot(nrow(points) == 13741L, all(points$qval < .05),
  all(is.finite(points$x)), all(is.finite(points$y)))
# Labels are a separate design input. Filter them after the novel targets exist.
label_path <- prep_input("paper", "output/qtl_circos_simplified_options/L_noncoding_highlights/gene_labels.tsv")
requested_labels <- read_tsv(label_path, show_col_types = FALSE) %>%
  semi_join(points %>% distinct(genes = gene_name), by = "genes")
stopifnot(nrow(requested_labels) == 26)
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
leaders <- anchors %>% select(label, sx, sy, label_theta) %>%
  left_join(gene_type_dots %>% group_by(label) %>%
              summarise(end_radius = min(dot_radius) - 0.85, .groups = "drop"),
            by = "label", relationship = "one-to-one") %>%
  mutate(ex = end_radius * cos(label_theta * pi / 180),
         ey = end_radius * sin(label_theta * pi / 180)) %>%
  select(label, sx, sy, ex, ey)


lead_points <- points
highlights <- label_points
label_dots <- gene_type_dots
axis_ticks <- axes
message("Calculated 13,741 novel-trait points, 26 labels and all circos coordinates.")
