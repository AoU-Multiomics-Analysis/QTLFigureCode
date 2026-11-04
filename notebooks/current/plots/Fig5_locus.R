# Pure locus and prediction plot methods from Fig5/00_setup.R.
suppressPackageStartupMessages(library(scales))
CCRE_FILL <- c(PLS = '#E15759', ELS = '#F4A261')
nucleotide_colors <- c(A = '#E64B35', C = '#00A087', G = '#F39B7F', T = '#4DBBD5')
gene_models_for_window <- function(chr, win_start, win_end, models = list(spans = gene_spans, exons = gene_exons)) {
  spans <- models$spans %>% filter(chrom == chr, start <= win_end, end >= win_start)
  list(spans = spans, exons = models$exons %>% filter(chrom == chr, gene %in% spans$gene, start <= win_end, end >= win_start))
}

cs_palette <- function(levels) {
    levels <- unique(levels[!is.na(levels)])
    setNames(rcartocolor::carto_pal(max(3, length(levels)), "Prism")[seq_along(levels)], levels)
}

signed_square_trans <- function() scales::trans_new("signed_square", transform = function(x) sign(x) *
    x^2, inverse = function(x) sign(x) * sqrt(abs(x)))

fmt_af <- function(af) dplyr::case_when(is.na(af) ~ "NA", af >= 0.01 ~ sprintf("%.1f%%", af * 100), af >=
    0.001 ~ sprintf("%.2f%%", af * 100), TRUE ~ sprintf("%.3f%%", af * 100))

make_freq_label <- function(max_af, subpop, rare) {
    pop <- ifelse(is.na(subpop), "", sprintf(" (%s)", toupper(subpop)))
    base <- sprintf("MAF %s%s", fmt_af(max_af), pop)
    ifelse(rare, paste0("rare · ", base), base)
}

normalize_variant_id <- function(x) str_replace_all(str_replace_all(as.character(x), "[:*]", "_"), "_+",
    "_")

make_variant_id <- function(chrom, pos, ref, alt) normalize_variant_id(paste(chrom, pos, toupper(ref),
    toupper(alt), sep = "_"))

build_variant_freq <- function(eqtl) {
    eqtl %>% transmute(variant_key = normalize_variant_id(variant), max_af = gvs_max_af, subpop = gvs_max_subpop,
        rare, log2_aFC) %>% distinct(variant_key, .keep_all = TRUE) %>% mutate(freq_label = make_freq_label(max_af,
        subpop, rare), afc_label = sprintf("log₂ aFC %+.1f", log2_aFC))
}

ccre_bands <- function(ccre, scale = 1e+06, alpha = 0.16, vars = NULL) {
    if (is.null(ccre) || nrow(ccre) == 0)
        return(NULL)
    if (!is.null(vars))
        ccre <- filter(ccre, map2_lgl(start, end, function(s, e) any(vars >= s & vars <= e)))
    if (nrow(ccre) == 0)
        return(NULL)
    ccre <- select(mutate(ccre, fillc = if_else(grepl("PLS", type), CCRE_FILL["PLS"], CCRE_FILL["ELS"])),
        start, end, fillc)
    edges <- pivot_longer(ccre, c(start, end), values_to = "x")
    list(geom_rect(data = ccre, aes(xmin = start/scale, xmax = end/scale), ymin = -Inf, ymax = Inf, fill = ccre$fillc,
        alpha = alpha, inherit.aes = FALSE), geom_segment(data = edges, aes(x = x/scale, xend = x/scale),
        y = -Inf, yend = Inf, colour = edges$fillc, linewidth = 0.3, inherit.aes = FALSE))
}

gene_model_track <- function(xlim_mb, chrom, scale = 1e+06, exon_h = 0.34, label_size = 3.4, models = NULL,
    expand = expansion(mult = 0), x_title = "Genomic position (Mb)", x_labeller = waiver(), right_strip = TRUE) {
    win <- xlim_mb * scale
    if (is.null(models))
        models <- gene_models_for_window(chrom, win[1], win[2])
    sp <- models$spans
    ex <- filter(models$exons, gene %in% sp$gene)
    lab <- mutate(sp, tss = if_else(strand == "+", start, end), arr = if_else(strand == "+", 1, -1) *
        diff(win) * 0.045, tss_in = tss >= win[1] & tss <= win[2], tss_x = pmin(pmax(tss, win[1]), win[2]))
    tss_lab <- filter(lab, tss_in)
    p <- ggplot() + geom_segment(data = sp, aes(x = pmax(start/scale, xlim_mb[1]), xend = pmin(end/scale,
        xlim_mb[2]), y = 0, yend = 0), colour = "grey45", linewidth = 0.5) + geom_rect(data = ex, aes(xmin = start/scale,
        xmax = end/scale, ymin = -exon_h, ymax = exon_h), fill = "grey35") + geom_segment(data = tss_lab,
        aes(x = tss/scale, xend = tss/scale, y = 0, yend = exon_h * 2.2), colour = "grey45", linewidth = 0.5) +
        geom_segment(data = tss_lab, aes(x = tss/scale, xend = (tss + arr)/scale, y = exon_h * 2.2, yend = exon_h *
            2.2), colour = "grey45", linewidth = 0.5, arrow = arrow(length = unit(3.5, "pt"), type = "closed")) +
        geom_text(data = lab, aes(x = tss_x/scale, y = 0, label = gene), size = label_size, fontface = "italic",
            hjust = 1, nudge_x = -diff(xlim_mb) * 0.012) + scale_x_continuous(limits = xlim_mb, expand = expand,
        labels = x_labeller) + scale_y_continuous(limits = c(-exon_h * 1.6, exon_h * 2.9)) + coord_cartesian(clip = "off") +
        labs(x = x_title, y = NULL) + theme(panel.background = element_blank(), panel.grid = element_blank(),
        panel.border = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 11), axis.title.x = element_text(size = 12), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.line.y = element_blank(), axis.title.y = element_blank(),
        legend.position = "none", plot.margin = margin(1, 8, 4, 8))
    if (right_strip)
        p <- p + facet_wrap(~"x", strip.position = "right") + theme(strip.text = element_text(size = 12,
            colour = NA), strip.background = element_blank(), panel.spacing = unit(0, "pt"))
    p
}

regional_window_bp <- function(position, causal, zoom_center, zoom_window_bp, flank_bp = 5000) {
    focal <- position[causal & !is.na(position)]
    if (!length(focal))
        focal <- position
    win <- range(focal, na.rm = TRUE) + c(-1, 1) * flank_bp
    win <- c(min(win[1], zoom_center - zoom_window_bp), max(win[2], zoom_center + zoom_window_bp))
    tested <- range(position, na.rm = TRUE)
    c(max(win[1], tested[1]), min(win[2], tested[2]))
}

inset_bottom_clearing_points <- function(x_npc, y_npc, requested, inset_left, inset_right, marker_npc = 0.022,
    max_bottom = 0.4) {
    under <- y_npc[is.finite(x_npc) & is.finite(y_npc) & x_npc >= inset_left & x_npc <= inset_right]
    bottom <- requested
    for (y in sort(under)) if (abs(y - bottom) < marker_npc)
        bottom <- y + marker_npc
    min(bottom, max_bottom)
}

plot_finemapping_with_inset <- function(finemapping_df, pip_thresh = 0.9, zoom_center, zoom_window_bp = 500,
    inset_left = 0.24, inset_bottom = 0.05, inset_right = 0.96, inset_top = 0.76, box_padding_bp = 1000,
    box_padding_pip = 0.05, cs_colors, panel_title = NULL, ccre = NULL, x_limits = NULL, chrom = NULL,
    regional_flank_bp = 5000, edge_pad_x = 0.015) {
    plot_df <- finemapping_df %>% filter(!is.na(position), !is.na(pip)) %>% mutate(causal = pip >= pip_thresh,
        cs_id = factor(cs_id), in_zoom = abs(position - zoom_center) <= zoom_window_bp)
    zoom_df <- plot_df %>% filter(in_zoom)
    iz_causal <- plot_df$position[plot_df$causal & plot_df$in_zoom]
    oz_causal <- plot_df$position[plot_df$causal & !plot_df$in_zoom]
    regional_x_limits <- if (is.null(x_limits))
        regional_window_bp(plot_df$position, plot_df$causal, zoom_center, zoom_window_bp, regional_flank_bp)/1e+06
    else x_limits
    regional_x_limits <- regional_x_limits + c(-1, 1) * diff(regional_x_limits) * edge_pad_x
    if (is.null(chrom)) {
        if (!"chrom" %in% names(plot_df))
            stop("plot_finemapping_with_inset(): pass chrom= (or supply a `chrom` ", "column) -- the gene-model track needs a chromosome to select genes.")
        chrom <- unique(plot_df$chrom)
        if (length(chrom) != 1L)
            stop("expected one chromosome, got: ", paste(chrom, collapse = ", "))
    }
    regional_y_limits <- c(0, 1.08)
    boxed <- zoom_df %>% filter(causal)
    if (nrow(boxed) == 0)
        boxed <- zoom_df
    box_start_mb <- max(min(boxed$position) - box_padding_bp, min(plot_df$position))/1e+06
    box_end_mb <- min(max(boxed$position) + box_padding_bp, max(plot_df$position))/1e+06
    box_center_mb <- mean(c(box_start_mb, box_end_mb))
    box_y_min <- max(regional_y_limits[1], min(boxed$pip) - box_padding_pip)
    box_y_max <- min(regional_y_limits[2], max(boxed$pip) + box_padding_pip)
    inset_bottom <- inset_bottom_clearing_points(x_npc = (plot_df$position/1e+06 - regional_x_limits[1])/diff(regional_x_limits),
        y_npc = (plot_df$pip - regional_y_limits[1])/diff(regional_y_limits), requested = inset_bottom,
        inset_left = inset_left, inset_right = inset_right)
    inset_top_data <- regional_y_limits[1] + inset_top * diff(regional_y_limits)
    regional_plot <- ggplot(plot_df, aes(position/1e+06, pip)) + ccre_bands(ccre, scale = 1e+06, vars = oz_causal) +
        geom_hline(yintercept = pip_thresh, linetype = "dashed", linewidth = 0.5) + annotate("segment",
        x = box_center_mb, xend = box_center_mb, y = box_y_min, yend = inset_top_data, color = "grey55",
        linewidth = 0.4) + annotate("rect", xmin = box_start_mb, xmax = box_end_mb, ymin = box_y_min,
        ymax = box_y_max, fill = "grey90", alpha = 0.18, color = "grey35", linewidth = 0.5) + geom_point(data = filter(plot_df,
        !causal), aes(fill = cs_id), shape = 21, size = 1.8, stroke = 0.35, color = "grey55", alpha = 0.5,
        show.legend = FALSE) + geom_point(data = filter(plot_df, causal), aes(fill = cs_id), shape = 23,
        size = 3.2, stroke = 0.65, color = "black", show.legend = FALSE) + geom_point(data = filter(plot_df,
        causal & rare & !in_zoom), shape = 5, size = 4.8, colour = "#8C6D1F", stroke = 0.9, show.legend = FALSE) +
        geom_point(data = tibble(type = plot_df$type[1]), aes(x = regional_x_limits[1], y = -0.5, shape = "rare (MAF < 1%)"),
            colour = "#8C6D1F", size = 3.2, stroke = 0.9, inherit.aes = FALSE) + scale_shape_manual(name = NULL,
        values = c(`rare (MAF < 1%)` = 5)) + facet_wrap(~type, ncol = 1, strip.position = "right") +
        scale_fill_manual(values = cs_colors, na.value = "grey80") + scale_x_continuous(limits = regional_x_limits,
        expand = expansion(mult = 0)) + scale_y_continuous(breaks = seq(0, 1, 0.25), limits = regional_y_limits,
        expand = expansion(mult = 0)) + labs(x = NULL, y = "PIP", title = panel_title) + theme(strip.text = element_text(size = 12),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.6), axis.line = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "top", legend.justification = "left",
        legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(1, "pt"), legend.text = element_text(size = 9),
        legend.key.size = unit(4, "mm"), legend.key = element_blank(), plot.title = element_text(size = 12.5,
            hjust = 0.5, face = "bold", margin = margin(b = 4)), plot.margin = margin(6, 8, 2, 8))
    zoom_plot <- ggplot(zoom_df, aes(position, pip)) + ccre_bands(ccre, scale = 1, vars = iz_causal) +
        geom_hline(yintercept = pip_thresh, linetype = "dashed", linewidth = 0.4) + geom_point(data = filter(zoom_df,
        !causal), aes(fill = cs_id), shape = 21, size = 1.7, stroke = 0.3, color = "grey55", alpha = 0.5,
        show.legend = FALSE) + geom_point(data = filter(zoom_df, causal), aes(fill = cs_id), shape = 23,
        size = 3, stroke = 0.6, color = "black", show.legend = FALSE) + geom_point(data = filter(zoom_df,
        causal & rare), shape = 5, size = 4.6, colour = "#8C6D1F", stroke = 0.85, show.legend = FALSE) +
        scale_fill_manual(values = cs_colors, na.value = "grey80") + scale_x_continuous(breaks = (function(b) b[seq(1,
        length(b), 2)])(pretty(zoom_center + c(-1, 1) * zoom_window_bp, 4)), limits = c(zoom_center -
        zoom_window_bp, zoom_center + zoom_window_bp), labels = function(x) format(x, big.mark = ",",
        scientific = FALSE, trim = TRUE), expand = expansion(mult = c(0.04, 0.04))) + scale_y_continuous(breaks = c(0,
        0.5, pip_thresh), limits = c(0, 1.12), expand = expansion(mult = c(0.02, 0.08))) + labs(x = NULL,
        y = NULL) + theme(panel.border = element_rect(fill = NA, color = "grey35", linewidth = 0.5),
        axis.line = element_blank(), axis.text = element_text(size = 8.5), legend.position = "none",
        plot.background = element_rect(fill = "white", color = NA))
    main <- regional_plot + inset_element(zoom_plot, left = inset_left, bottom = inset_bottom, right = inset_right,
        top = inset_top, align_to = "panel", on_top = TRUE, ignore_tag = TRUE)
    gtrack <- gene_model_track(regional_x_limits, chrom = chrom, expand = expansion(mult = 0))
    main/gtrack + plot_layout(heights = c(1, 0.17))
}

prepare_promoterAI_locus <- function(promoter_data, qtl_data, gene_symbol, transcript, focal_position,
    window_bp = 100) {
    qtl_lookup <- qtl_data %>% transmute(variant_key = normalize_variant_id(variant), pip, causal, cs_id) %>%
        distinct(variant_key, .keep_all = TRUE)
    promoter_data %>% filter(gene == gene_symbol, transcript_id == transcript) %>% mutate(relative_position = pos -
        focal_position, variant = make_variant_id(chrom, pos, ref, alt), variant_key = normalize_variant_id(variant)) %>%
        filter(abs(relative_position) < window_bp) %>% left_join(qtl_lookup, by = "variant_key") %>%
        mutate(causal = replace_na(causal, FALSE))
}

build_causal_promoterAI_facets <- function(promoter_data, cs_palette_vec, flank_bp = 20, freq_lookup = NULL) {
    plot_data <- promoter_data %>% mutate(across(c(chrom, ref, alt, cs_id), as.character), pos = as.numeric(pos),
        causal = coalesce(as.logical(causal), FALSE), variant_label = paste0(chrom, ":", format(pos,
            scientific = FALSE, trim = TRUE), " ", ref, ">", alt))
    causal_variants <- plot_data %>% filter(causal) %>% distinct(chrom, pos, ref, alt, cs_id, variant_label) %>%
        arrange(chrom, pos, ref, alt) %>% mutate(finemap_color = unname(cs_palette_vec[cs_id]), variant_key = normalize_variant_id(make_variant_id(chrom,
        pos, ref, alt)))
    if (!is.null(freq_lookup))
        causal_variants <- causal_variants %>% left_join(freq_lookup %>% select(variant_key, freq_label,
            afc_label), by = "variant_key") %>% mutate(strip_label = if_else(is.na(freq_label), variant_label,
            paste0(variant_label, "\n", freq_label, "\n", afc_label)))
    else causal_variants <- causal_variants %>% mutate(strip_label = variant_label)
    facet_information <- causal_variants %>% transmute(chrom, facet_id = paste(chrom, format(pos, scientific = FALSE,
        trim = TRUE), ref, alt, sep = "_"), focal_position = pos, focal_ref = ref, focal_alt = alt, window_start = pos -
        flank_bp, window_end = pos + flank_bp, cs_id, finemap_color, facet_strip_label = strip_label) %>%
        arrange(chrom, focal_position) %>% mutate(facet_strip_label = factor(facet_strip_label, levels = unique(facet_strip_label)))
    facet_plot_data <- plot_data %>% inner_join(facet_information, by = "chrom", relationship = "many-to-many") %>%
        filter(pos >= window_start, pos <= window_end) %>% mutate(relative_position = pos - focal_position,
        is_facet_causal = pos == focal_position & ref == focal_ref & alt == focal_alt, facet_strip_label = factor(facet_strip_label,
            levels = levels(facet_information$facet_strip_label)))
    attr(facet_plot_data, "facet_information") <- facet_information
    facet_plot_data
}

plot_causal_promoterAI_facets <- function(promoter_data, cs_palette_vec, flank_bp = 20, effect_threshold = 0.4,
    causal_point_size = 4, facet_diamond_size = 3, facet_diamond_stroke = 0.8, motif_annotations = NULL,
    freq_lookup = NULL, nucleotide_cols = nucleotide_colors) {
    fpd <- build_causal_promoterAI_facets(promoter_data, cs_palette_vec, flank_bp, freq_lookup)
    facet_info <- attr(fpd, "facet_information")
    causal_pts <- fpd %>% filter(is_facet_causal) %>% distinct(facet_id, facet_strip_label, chrom, pos,
        ref, alt, promoterAI, relative_position, .keep_all = TRUE)
    facet_symbols <- facet_info %>% transmute(facet_id, facet_strip_label, relative_position = 0, finemap_color)
    y_limit <- max(1, ceiling(max(abs(fpd$promoterAI), na.rm = TRUE) * 10)/10)
    y_breaks <- sort(unique(c(-y_limit, -0.75, -0.5, 0, 0.5, 0.75, y_limit)))
    y_breaks <- y_breaks[y_breaks >= -y_limit & y_breaks <= y_limit]
    motif_col <- "#3E86C0"
    motif_layers <- if (!is.null(motif_annotations)) {
        ann <- motif_annotations %>% left_join(facet_info %>% select(facet_id, facet_strip_label), by = "facet_id") %>%
            mutate(facet_strip_label = factor(facet_strip_label, levels = levels(fpd$facet_strip_label)))
        list(geom_rect(data = ann, aes(xmin = xmin, xmax = xmax), ymin = -Inf, ymax = y_limit, fill = motif_col,
            alpha = 0.22, inherit.aes = FALSE), geom_text(data = ann, aes(x = (xmin + xmax)/2, label = paste0(tf,
            "\nmotif")), y = y_limit * 0.82, inherit.aes = FALSE, fontface = "bold.italic", size = 3.2,
            colour = "grey15", lineheight = 0.85))
    }
    else NULL
    ggplot(fpd, aes(relative_position, promoterAI)) + motif_layers + annotate("segment", x = -flank_bp *
        0.88, xend = flank_bp * 0.88, y = y_limit * 1.12, yend = y_limit * 1.12, colour = "grey72", linewidth = 0.6,
        lineend = "round") + geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.45) + geom_hline(yintercept = c(-effect_threshold,
        effect_threshold), linetype = "dotted", color = "grey68", linewidth = 0.4) + geom_segment(aes(xend = relative_position,
        y = 0, yend = promoterAI), color = "grey72", linewidth = 0.35, alpha = 0.75) + geom_point(data = filter(fpd,
        !is_facet_causal, abs(promoterAI) < effect_threshold), shape = 16, size = 1.5, color = "grey76",
        alpha = 0.55, show.legend = FALSE) + geom_point(data = filter(fpd, !is_facet_causal, abs(promoterAI) >=
        effect_threshold), aes(color = alt), shape = 16, size = 2.3, alpha = 0.9) + geom_point(data = causal_pts,
        aes(fill = alt), shape = 21, size = causal_point_size, stroke = 0.7, color = "black", show.legend = FALSE) +
        geom_point(data = facet_symbols, aes(relative_position, y = Inf, fill = I(finemap_color)), inherit.aes = FALSE,
            shape = 23, size = facet_diamond_size, stroke = facet_diamond_stroke, color = "black", show.legend = FALSE) +
        facet_grid(cols = vars(facet_strip_label), scales = "free_x", space = "free_x") + scale_color_manual(values = nucleotide_cols,
        breaks = c("A", "C", "T", "G"), drop = FALSE, name = "Alternate allele") + scale_fill_manual(values = nucleotide_cols,
        breaks = c("A", "C", "T", "G"), guide = "none") + scale_x_continuous(breaks = function(x) pretty(x,
        5), expand = expansion(mult = c(0.025, 0.025))) + scale_y_continuous(trans = signed_square_trans(),
        limits = c(-y_limit, y_limit), oob = scales::oob_keep, breaks = y_breaks, labels = label_number(accuracy = 0.01,
            trim = TRUE), expand = expansion(mult = c(0.05, 0.15))) + coord_cartesian(clip = "off") +
        labs(x = "Position relative to causal variant (bp)", y = "PromoterAI effect") + guides(color = guide_legend(title.position = "top",
        ncol = 2, byrow = TRUE, override.aes = list(size = 2.5, alpha = 1))) + theme(strip.text.x = element_text(size = 9,
        face = "bold", margin = margin(2, 2, 15, 2)), panel.spacing.x = unit(1.2, "lines"), axis.text = element_text(size = 9),
        axis.title = element_text(size = 11), legend.position = "inside", legend.position.inside = c(0.02,
            0.04), legend.justification = c(0, 0), legend.title = element_text(size = 7), legend.text = element_text(size = 7),
        legend.key.size = unit(3.5, "mm"), legend.background = element_rect(fill = alpha("white", 0.85),
            color = NA), plot.margin = margin(8, 15, 5, 5))
}

.acvrl1_qtl <- function(eqtl) eqtl %>% filter(gene_name == "ACVRL1") %>% mutate(causal = pip > 0.9, cs_id = as.character(cs_id))

# The locus strip identifies the molecular assay, as in the original panel D.
fig5_plot_acvrl1_locus <- function(qtl, ccre, focal_position = 51912316) {
  plot_finemapping_with_inset(qtl %>% mutate(type = "Expression"),
    pip_thresh = 0.9, zoom_center = focal_position, zoom_window_bp = 500,
    cs_colors = cs_palette(qtl$cs_id), chrom = "chr12", panel_title = NULL, ccre = ccre)
}

# Select GENCODE gene and exon rows for this genomic window. Merge overlapping
# and adjacent exons within each gene, as in common/gene_models.R.
s31_read_gene_models <- function(path, chromosome, window_start, window_end) {
  connection <- if (grepl("[.]gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(connection))
  chunks <- list()
  repeat {
    lines <- readLines(connection, n = 50000L, warn = FALSE)
    if (!length(lines)) break
    lines <- lines[str_starts(lines, fixed(paste0(chromosome, "\t"))) &
                     str_detect(lines, "^[^\t]+\t[^\t]+\t(gene|exon)\t")]
    if (!length(lines)) next
    rows <- read_tsv(I(paste(lines, collapse = "\n")),
      col_names = c("chrom", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"),
      col_types = cols(.default = col_character()), quote = "", progress = FALSE) %>%
      transmute(gene = str_match(attributes, 'gene_name "([^"]+)"')[, 2], chrom,
        start = as.integer(start), end = as.integer(end), strand, feature,
        gene_type = str_match(attributes, 'gene_type "([^"]+)"')[, 2]) %>%
      filter(start <= window_end, end >= window_start, gene_type == "protein_coding", !is.na(gene))
    chunks[[length(chunks) + 1L]] <- rows
  }
  rows <- bind_rows(chunks)
  spans <- rows %>% filter(feature == "gene") %>%
    arrange(gene, desc(end - start)) %>% distinct(gene, .keep_all = TRUE) %>%
    select(gene, chrom, start, end, strand) %>% arrange(chrom, start)
  exons <- rows %>% filter(feature == "exon", gene %in% spans$gene) %>%
    group_by(gene, chrom, strand) %>% arrange(start, end, .by_group = TRUE) %>%
    mutate(block = cumsum(start > lag(as.double(cummax(end)), default = -Inf) + 1)) %>%
    group_by(gene, chrom, strand, block) %>%
    summarise(start = min(start), end = max(end), .groups = "drop") %>%
    select(gene, chrom, start, end, strand) %>% arrange(chrom, start)
  list(spans = spans, exons = exons)
}
