## ============================================================================
## Figure 5 panels D/E — ACVRL1 fine-mapping locus + PromoterAI base-resolution
## Sourced by PlotFigure5.R (expects `eQTL_susie` and make_signed_square_trans()
## in scope).  Produces: acvrl1_finemap_plot, acvrl1_promoterAI_plot
## ============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(data.table)
  library(patchwork); library(ggrepel); library(rcartocolor)
})

ACVRL1_FOCAL <- 51912316   # ACVRL1 promoter variant (chr12, hg38)

## ---- data -----------------------------------------------------------------
ACVRL1_QTL_data <- eQTL_susie %>%
  filter(gene_name == "ACVRL1") %>%
  mutate(causal = pip > 0.9, cs_id = as.character(cs_id))

promoterAI_data <- fread(file.path(DIR_FIG, "acvrl1_promoterAI.tsv"))
acvrl1_ccre     <- fread(file.path(DIR_FIG, "acvrl1_ccre_region.tsv"))  # covers the full D locus

## ENCODE cCRE shading: PLS = promoter-like (red), (d/p)ELS = enhancer-like (orange)
CCRE_FILL <- c(PLS = "#E15759", ELS = "#F4A261")
ccre_bands <- function(ccre, scale = 1e6, alpha = 0.16, vars = NULL) {
  if (is.null(ccre) || nrow(ccre) == 0) return(NULL)
  ## keep only cCREs overlapping a fine-mapped variant position
  if (!is.null(vars))
    ccre <- ccre |> filter(map2_lgl(start, end, \(s, e) any(vars >= s & vars <= e)))
  if (nrow(ccre) == 0) return(NULL)
  ## drop `type` (would collide with the plots' facet_wrap(~type)); keep only geometry + fill
  ccre <- ccre |>
    mutate(fillc = if_else(grepl("PLS", type), CCRE_FILL["PLS"], CCRE_FILL["ELS"])) |>
    select(start, end, fillc)
  ## light fill + crisp vertical edges only (no top/bottom rule) so elements stay distinct
  edges <- ccre |> pivot_longer(c(start, end), values_to = "x")
  list(
    geom_rect(data = ccre, aes(xmin = start / scale, xmax = end / scale),
              ymin = -Inf, ymax = Inf, fill = ccre$fillc, alpha = alpha, inherit.aes = FALSE),
    geom_segment(data = edges, aes(x = x / scale, xend = x / scale), y = -Inf, yend = Inf,
                 colour = edges$fillc, linewidth = 0.3, inherit.aes = FALSE)
  )
}

## ---- gene-model track (aligns under a fine-map panel via shared x limits) ----
gene_exons_all <- fread(file.path(DIR_FIG, "gene_exons.tsv"),
                        col.names = c("gene", "chrom", "start", "end", "strand"))
gene_spans_all <- fread(file.path(DIR_FIG, "gene_spans.tsv"),
                        col.names = c("gene", "chrom", "start", "end", "strand"))

gene_model_track <- function(xlim_mb, scale = 1e6, exon_h = 0.34, label_size = 3.4,
                             exons = gene_exons_all, spans = gene_spans_all,
                             expand = expansion(mult = 0), x_title = "Genomic position (Mb)",
                             x_labeller = waiver(), right_strip = TRUE) {
  win <- xlim_mb * scale
  sp  <- spans |> filter(start < win[2], end > win[1])
  ex  <- exons |> filter(gene %in% sp$gene, start < win[2], end > win[1])
  lab <- sp |> mutate(vmid = (pmax(start, win[1]) + pmin(end, win[2])) / 2,
                      tss  = if_else(strand == "+", start, end),
                      arr  = if_else(strand == "+", 1, -1) * diff(win) * 0.045)
  p <- ggplot() +
    ## gene-body/intron line, clamped to the window so genes running off-panel aren't dropped
    geom_segment(data = sp, aes(x = pmax(start / scale, xlim_mb[1]),
                                xend = pmin(end / scale, xlim_mb[2]), y = 0, yend = 0),
                 colour = "grey45", linewidth = 0.5) +
    geom_rect(data = ex, aes(xmin = start / scale, xmax = end / scale,
                             ymin = -exon_h, ymax = exon_h),
              fill = "grey35") +
    geom_segment(data = lab, aes(x = tss / scale, xend = tss / scale,
                                 y = 0, yend = exon_h * 2.2), colour = "grey45", linewidth = 0.5) +
    geom_segment(data = lab, aes(x = tss / scale, xend = (tss + arr) / scale,
                                 y = exon_h * 2.2, yend = exon_h * 2.2), colour = "grey45",
                 linewidth = 0.5, arrow = arrow(length = unit(3.5, "pt"), type = "closed")) +
    geom_text(data = lab, aes(x = tss / scale, y = 0, label = gene),   # left of TSS, in line
              size = label_size, fontface = "italic", hjust = 1,
              nudge_x = -diff(xlim_mb) * 0.012) +
    scale_x_continuous(limits = xlim_mb, expand = expand, labels = x_labeller) +
    scale_y_continuous(limits = c(-exon_h * 1.6, exon_h * 2.9)) +
    coord_cartesian(clip = "off") +
    labs(x = x_title, y = NULL) +
    theme(panel.background = element_blank(), panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line.x = element_blank(), axis.ticks.x = element_blank(),  # gene line is the only rule
          axis.text.x = element_text(size = 11), axis.title.x = element_text(size = 12),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.line.y = element_blank(), axis.title.y = element_blank(),
          legend.position = "none", plot.margin = margin(1, 8, 4, 8))
  ## matching (invisible) right strip keeps the panel width aligned with the faceted plot above
  if (right_strip)
    p <- p + facet_wrap(~"x", strip.position = "right") +
      theme(strip.text = element_text(size = 12, colour = NA),
            strip.background = element_blank(), panel.spacing = unit(0, "pt"))
  p
}

## shared credible-set palette (Prism), used by both D and E
acvrl1_cs_levels <- ACVRL1_QTL_data %>% filter(!is.na(cs_id)) %>% pull(cs_id) %>% unique()
shared_cs_palette <- setNames(
  carto_pal(max(3, length(acvrl1_cs_levels)), "Prism")[seq_along(acvrl1_cs_levels)],
  acvrl1_cs_levels)

## ===========================================================================
## Panel D — regional PIP plot with a zoom inset
## ===========================================================================
plot_finemapping_with_inset <- function(finemapping_df, pip_thresh = 0.9,
    zoom_center, zoom_window_bp = 500,
    inset_left = 0.24, inset_bottom = 0.05, inset_right = 0.96, inset_top = 0.76,
    box_padding_bp = 1000, box_padding_pip = 0.05, cs_colors, panel_tag = NULL,
    panel_title = NULL, ccre = NULL, x_limits = NULL) {

  plot_df <- finemapping_df %>%
    filter(!is.na(position), !is.na(pip)) %>%
    mutate(causal = pip >= pip_thresh, cs_id = factor(cs_id),
           in_zoom = abs(position - zoom_center) <= zoom_window_bp)
  zoom_df <- plot_df %>% filter(in_zoom)
  iz_causal <- plot_df$position[plot_df$causal &  plot_df$in_zoom]   # inset gets these cCREs/rings
  oz_causal <- plot_df$position[plot_df$causal & !plot_df$in_zoom]   # regional flags these

  regional_x_range <- range(plot_df$position / 1e6, na.rm = TRUE)
  regional_x_pad   <- max(diff(regional_x_range) * 0.02, 0.001)
  regional_x_limits<- if (is.null(x_limits)) regional_x_range + c(-1, 1) * regional_x_pad else x_limits
  regional_y_limits<- c(0, 1.08)

  boxed <- zoom_df %>% filter(causal); if (nrow(boxed) == 0) boxed <- zoom_df
  box_start_mb <- max(min(boxed$position) - box_padding_bp, min(plot_df$position)) / 1e6
  box_end_mb   <- min(max(boxed$position) + box_padding_bp, max(plot_df$position)) / 1e6
  box_center_mb<- mean(c(box_start_mb, box_end_mb))
  box_y_min <- max(regional_y_limits[1], min(boxed$pip) - box_padding_pip)
  box_y_max <- min(regional_y_limits[2], max(boxed$pip) + box_padding_pip)
  inset_top_data <- regional_y_limits[1] + inset_top * diff(regional_y_limits)

  regional_plot <- ggplot(plot_df, aes(position / 1e6, pip)) +
    ## regulatory-element highlight for causal variants that fall outside the inset
    ccre_bands(ccre, scale = 1e6, vars = oz_causal) +
    geom_hline(yintercept = pip_thresh, linetype = "dashed", linewidth = 0.5) +
    annotate("segment", x = box_center_mb, xend = box_center_mb,
             y = box_y_min, yend = inset_top_data, color = "grey55", linewidth = 0.4) +
    annotate("rect", xmin = box_start_mb, xmax = box_end_mb, ymin = box_y_min, ymax = box_y_max,
             fill = "grey90", alpha = 0.18, color = "grey35", linewidth = 0.5) +
    geom_point(data = filter(plot_df, !causal), aes(fill = cs_id), shape = 21,
               size = 1.8, stroke = 0.35, color = "grey55", alpha = 0.5, show.legend = FALSE) +
    geom_point(data = filter(plot_df, causal), aes(fill = cs_id), shape = 23,
               size = 3.2, stroke = 0.65, color = "black", show.legend = FALSE) +
    ## rare ring for out-of-zoom causal variants (in-zoom ones are ringed in the inset)
    geom_point(data = filter(plot_df, causal & rare & !in_zoom), shape = 5, size = 4.8,
               colour = "#8C6D1F", stroke = 0.9, show.legend = FALSE) +
    ## rare-variant key rendered as a top-left legend (an off-panel dummy point drives it)
    geom_point(data = tibble(type = plot_df$type[1]),
               aes(x = regional_x_limits[1], y = -0.5, shape = "rare (MAF < 1%)"),
               colour = "#8C6D1F", size = 3.2, stroke = 0.9, inherit.aes = FALSE) +
    scale_shape_manual(name = NULL, values = c("rare (MAF < 1%)" = 5)) +
    facet_wrap(~type, ncol = 1, strip.position = "right") +
    scale_fill_manual(values = cs_colors, na.value = "grey80") +
    scale_x_continuous(limits = regional_x_limits, expand = expansion(mult = 0)) +
    scale_y_continuous(breaks = seq(0, 1, 0.25), limits = regional_y_limits,
                       expand = expansion(mult = 0)) +
    labs(x = NULL, y = "PIP", tag = panel_tag, title = panel_title) +
    theme(strip.text = element_text(size = 12),
          panel.border = element_rect(fill = NA, color = "black", linewidth = 0.6),
          axis.line = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),  # gene track owns the axis
          legend.position = "top", legend.justification = "left",
          legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(1, "pt"),
          legend.text = element_text(size = 9), legend.key.size = unit(4, "mm"),
          legend.key = element_blank(),
          plot.title = element_text(size = 12.5, hjust = 0.5, face = "bold",
                                    margin = margin(b = 4)),
          plot.tag.position = c(0.02, 1.02), plot.margin = margin(6, 8, 2, 8))

  zoom_plot <- ggplot(zoom_df, aes(position, pip)) +
    ccre_bands(ccre, scale = 1, vars = iz_causal) +
    geom_hline(yintercept = pip_thresh, linetype = "dashed", linewidth = 0.4) +
    geom_point(data = filter(zoom_df, !causal), aes(fill = cs_id), shape = 21,
               size = 1.7, stroke = 0.3, color = "grey55", alpha = 0.5, show.legend = FALSE) +
    geom_point(data = filter(zoom_df, causal), aes(fill = cs_id), shape = 23,
               size = 3, stroke = 0.6, color = "black", show.legend = FALSE) +
    geom_point(data = filter(zoom_df, causal & rare), shape = 5, size = 4.6,
               colour = "#8C6D1F", stroke = 0.85, show.legend = FALSE) +
    scale_fill_manual(values = cs_colors, na.value = "grey80") +
    scale_x_continuous(                                    # every other tick, to de-clutter
      breaks = (\(b) b[seq(1, length(b), 2)])(pretty(zoom_center + c(-1, 1) * zoom_window_bp, 4)),
      limits = c(zoom_center - zoom_window_bp, zoom_center + zoom_window_bp),
      labels = \(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE),
      expand = expansion(mult = c(0.04, 0.04))) +
    scale_y_continuous(breaks = c(0, 0.5, pip_thresh), limits = c(0, 1.12),
                       expand = expansion(mult = c(0.02, 0.08))) +
    labs(x = NULL, y = NULL) +          # inset: drop axis titles
    theme(panel.border = element_rect(fill = NA, color = "grey35", linewidth = 0.5),
          axis.line = element_blank(), axis.text = element_text(size = 8.5),
          legend.position = "none",
          plot.background = element_rect(fill = "white", color = NA))

  main <- regional_plot +
    inset_element(zoom_plot, left = inset_left, bottom = inset_bottom,
                  right = inset_right, top = inset_top, align_to = "panel", on_top = TRUE)
  gtrack <- gene_model_track(regional_x_limits, expand = expansion(mult = 0))
  main / gtrack + plot_layout(heights = c(1, 0.17))
}

acvrl1_finemap_plot <- plot_finemapping_with_inset(
  ACVRL1_QTL_data %>% mutate(type = "Expression"),
  pip_thresh = 0.9, zoom_center = ACVRL1_FOCAL, zoom_window_bp = 500,
  cs_colors = shared_cs_palette, panel_tag = "D", panel_title = "ACVRL1 eQTL fine-mapping",
  ccre = acvrl1_ccre)

## ===========================================================================
## Panel E — PromoterAI base-resolution effect, one facet per causal variant
## ===========================================================================
normalize_variant_id <- function(x) str_replace_all(str_replace_all(as.character(x), "[:*]", "_"), "_+", "_")
make_variant_id <- function(chrom, pos, ref, alt)
  normalize_variant_id(paste(chrom, pos, toupper(ref), toupper(alt), sep = "_"))

## PromoterAI for one transcript around a focal position, joined to QTL
## (now carries cs_id + causal through, which the facet builder needs)
prepare_promoterAI_locus <- function(promoter_data, qtl_data, gene_symbol, transcript,
                                     focal_position, window_bp = 100) {
  qtl_lookup <- qtl_data %>%
    transmute(variant_key = normalize_variant_id(variant), pip, causal, cs_id) %>%
    distinct(variant_key, .keep_all = TRUE)
  promoter_data %>%
    filter(gene == gene_symbol, transcript_id == transcript) %>%
    mutate(relative_position = pos - focal_position,
           variant     = make_variant_id(chrom, pos, ref, alt),
           variant_key = normalize_variant_id(variant)) %>%
    filter(abs(relative_position) < window_bp) %>%
    left_join(qtl_lookup, by = "variant_key") %>%
    mutate(causal = replace_na(causal, FALSE))
}

## one independent window per causal variant
build_causal_promoterAI_facets <- function(promoter_data, cs_palette, flank_bp = 20,
                                           freq_lookup = NULL) {
  plot_data <- promoter_data %>%
    mutate(across(c(chrom, ref, alt, cs_id), as.character), pos = as.numeric(pos),
           causal = coalesce(as.logical(causal), FALSE),
           variant_label = paste0(chrom, ":", format(pos, scientific = FALSE, trim = TRUE),
                                   " ", ref, ">", alt))
  causal_variants <- plot_data %>% filter(causal) %>%
    distinct(chrom, pos, ref, alt, cs_id, variant_label) %>% arrange(chrom, pos, ref, alt) %>%
    mutate(finemap_color = unname(cs_palette[cs_id]),
           variant_key = normalize_variant_id(make_variant_id(chrom, pos, ref, alt)))
  ## append rarity / ancestry as a second strip line when a lookup is supplied
  if (!is.null(freq_lookup))
    causal_variants <- causal_variants %>%
      left_join(freq_lookup %>% select(variant_key, freq_label, afc_label), by = "variant_key") %>%
      mutate(strip_label = if_else(is.na(freq_label), variant_label,
                                   paste0(variant_label, "\n", freq_label, "\n", afc_label)))
  else
    causal_variants <- causal_variants %>% mutate(strip_label = variant_label)
  facet_information <- causal_variants %>%
    transmute(chrom,
              facet_id = paste(chrom, format(pos, scientific = FALSE, trim = TRUE), ref, alt, sep = "_"),
              focal_position = pos, focal_ref = ref, focal_alt = alt,
              window_start = pos - flank_bp, window_end = pos + flank_bp,
              cs_id, finemap_color,
              facet_strip_label = strip_label) %>%
    arrange(chrom, focal_position) %>%
    mutate(facet_strip_label = factor(facet_strip_label, levels = unique(facet_strip_label)))
  facet_plot_data <- plot_data %>%
    inner_join(facet_information, by = "chrom", relationship = "many-to-many") %>%
    filter(pos >= window_start, pos <= window_end) %>%
    mutate(relative_position = pos - focal_position,
           is_facet_causal = pos == focal_position & ref == focal_ref & alt == focal_alt,
           facet_strip_label = factor(facet_strip_label, levels = levels(facet_information$facet_strip_label)))
  attr(facet_plot_data, "facet_information") <- facet_information
  facet_plot_data
}

plot_causal_promoterAI_facets <- function(promoter_data, cs_palette, flank_bp = 20,
    effect_threshold = 0.4, causal_point_size = 4, facet_diamond_size = 3,
    facet_diamond_stroke = 0.8, motif_annotations = NULL, freq_lookup = NULL,
    nucleotide_colors = c(A = "#E64B35", C = "#00A087", G = "#F39B7F", T = "#4DBBD5")) {

  fpd <- build_causal_promoterAI_facets(promoter_data, cs_palette, flank_bp, freq_lookup)
  facet_info <- attr(fpd, "facet_information")
  causal_pts <- fpd %>% filter(is_facet_causal) %>%
    distinct(facet_id, facet_strip_label, chrom, pos, ref, alt, promoterAI, relative_position, .keep_all = TRUE)
  facet_symbols <- facet_info %>% transmute(facet_id, facet_strip_label, relative_position = 0, finemap_color)

  y_limit <- max(1, ceiling(max(abs(fpd$promoterAI), na.rm = TRUE) * 10) / 10)
  y_breaks <- sort(unique(c(-y_limit, -0.75, -0.5, 0, 0.5, 0.75, y_limit)))
  y_breaks <- y_breaks[y_breaks >= -y_limit & y_breaks <= y_limit]

  ## TF-motif highlight + label, drawn only in the annotated facet(s).
  ## Light-blue box (distinct from the orange/red cCRE annotations) with a crisp edge.
  motif_col <- "#3E86C0"
  motif_layers <- if (!is.null(motif_annotations)) {
    ## key by stable facet_id, then resolve the (possibly multi-line) strip label
    ann <- motif_annotations %>%
      left_join(facet_info %>% select(facet_id, facet_strip_label), by = "facet_id") %>%
      mutate(facet_strip_label = factor(facet_strip_label, levels = levels(fpd$facet_strip_label)))
    list(
      ## box top capped at y_limit so it doesn't run into the grey facet-header bar
      geom_rect(data = ann, aes(xmin = xmin, xmax = xmax), ymin = -Inf, ymax = y_limit,
                fill = motif_col, alpha = 0.22, inherit.aes = FALSE),
      geom_text(data = ann, aes(x = (xmin + xmax) / 2, label = paste0(tf, "\nmotif")),
                y = y_limit * 0.82, inherit.aes = FALSE, fontface = "bold.italic",
                size = 3.2, colour = "grey15", lineheight = 0.85)
    )
  } else NULL

  ggplot(fpd, aes(relative_position, promoterAI)) +
    motif_layers +
    ## grey bar beneath each facet label, sitting just below the diamond (panel B style)
    annotate("segment", x = -flank_bp * 0.88, xend = flank_bp * 0.88,
             y = y_limit * 1.12, yend = y_limit * 1.12,
             colour = "grey72", linewidth = 0.6, lineend = "round") +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.45) +
    geom_hline(yintercept = c(-effect_threshold, effect_threshold),
               linetype = "dotted", color = "grey68", linewidth = 0.4) +
    geom_segment(aes(xend = relative_position, y = 0, yend = promoterAI),
                 color = "grey72", linewidth = 0.35, alpha = 0.75) +
    geom_point(data = filter(fpd, !is_facet_causal, abs(promoterAI) < effect_threshold),
               shape = 16, size = 1.5, color = "grey76", alpha = 0.55, show.legend = FALSE) +
    geom_point(data = filter(fpd, !is_facet_causal, abs(promoterAI) >= effect_threshold),
               aes(color = alt), shape = 16, size = 2.3, alpha = 0.9) +
    geom_point(data = causal_pts, aes(fill = alt), shape = 21, size = causal_point_size,
               stroke = 0.7, color = "black", show.legend = FALSE) +
    geom_point(data = facet_symbols, aes(relative_position, y = Inf, fill = I(finemap_color)),
               inherit.aes = FALSE, shape = 23, size = facet_diamond_size,
               stroke = facet_diamond_stroke, color = "black", show.legend = FALSE) +
    facet_grid(cols = vars(facet_strip_label), scales = "free_x", space = "free_x") +
    scale_color_manual(values = nucleotide_colors, breaks = c("A","C","T","G"),
                       drop = FALSE, name = "Alternate allele") +
    scale_fill_manual(values = nucleotide_colors, breaks = c("A","C","T","G"), guide = "none") +
    scale_x_continuous(breaks = \(x) pretty(x, 5), expand = expansion(mult = c(0.025, 0.025))) +
    scale_y_continuous(trans = make_signed_square_trans(), limits = c(-y_limit, y_limit),
                       oob = scales::oob_keep,   # keep the header bar (drawn just above y_limit)
                       breaks = y_breaks, labels = label_number(accuracy = 0.01, trim = TRUE),
                       expand = expansion(mult = c(0.05, 0.15))) +
    coord_cartesian(clip = "off") +
    labs(x = "Position relative to causal variant (bp)", y = "PromoterAI effect") +
    guides(color = guide_legend(title.position = "top", ncol = 2, byrow = TRUE,
                                override.aes = list(size = 2.5, alpha = 1))) +
    theme(strip.text.x = element_text(size = 9, face = "bold",
                                      margin = margin(2, 2, 15, 2)),
          panel.spacing.x = unit(1.2, "lines"),
          axis.text = element_text(size = 9), axis.title = element_text(size = 11),
          legend.position = "inside", legend.position.inside = c(0.02, 0.04),
          legend.justification = c(0, 0), legend.title = element_text(size = 7),
          legend.text = element_text(size = 7), legend.key.size = unit(3.5, "mm"),
          legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
          plot.margin = margin(8, 15, 5, 5))
}

ACVRL1_promoterAI <- prepare_promoterAI_locus(
  promoterAI_data, ACVRL1_QTL_data,
  gene_symbol = "ACVRL1", transcript = "ENST00000419526.6",
  focal_position = ACVRL1_FOCAL, window_bp = 100)

## focal green-diamond variant (chr12:51912316 C>T) sits in an SP1 GC-box:
## PromoterAI -3..+5 = GCCCCGCCC, whose rev-comp GGGCGGGGC is the SP1 consensus
acvrl1_promoterAI_plot <- plot_causal_promoterAI_facets(
  ACVRL1_promoterAI, cs_palette = shared_cs_palette, flank_bp = 20, effect_threshold = 0.5,
  freq_lookup = variant_freq,
  motif_annotations = tibble(facet_id = "chr12_51912316_C_T",
                             xmin = -3.5, xmax = 5.5, tf = "SP1"))
