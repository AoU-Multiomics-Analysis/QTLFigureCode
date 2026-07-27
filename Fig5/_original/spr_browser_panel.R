## ============================================================================
## SPR Borzoi multi-track browser panel  (Figure 5, second worked example)
##
## Self-contained: source this file, then call make_spr_browser() to get a
## cowplot object you can embed with patchwork::wrap_elements().
##
## Inputs (absolute paths): Borzoi-predicted ref/alt bigWigs and the CTCF
## gradient x input attribution table, generated in ~/borzoi_run.
## ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
  library(cowplot)
  library(ggseqlogo)
})

## ---- config (absolute so it works from any working directory) --------------
SPR_BW_DIR   <- "/Users/evinmpadhi/borzoi_run/exports/SPR_bigwig"
SPR_ATTR     <- "/Users/evinmpadhi/borzoi_run/results/attributions/spr_attr.tsv"
SPR_VARIANT_POS <- 72887998   # SPR variant, chr2 (hg38)

MARK_ORDER <- c("CTCF","cohesin","ELF1","accessibility","DNase","ATAC","EP300","H3K27ac",
                "H3K4me1","H3K4me3","CAGE","RNA","H3K36me3","EZH2","H3K27me3","H3K9me3")

## per-mark reference-track colours -- cool NASA "Visions of the Future" tones,
## kept red-free so the red alternate line always contrasts.
MARK_COLORS <- c(
  CTCF = "#253473", cohesin = "#3CC8C0", ELF1 = "#61619B", accessibility = "#67A850",
  DNase = "#67A850", ATAC = "#41AE76",
  EP300 = "#8CCCE4", H3K27ac = "#4982B6", H3K4me1 = "#4E9678", H3K4me3 = "#303251",
  CAGE = "#A0B8C3", RNA = "#2D332D", H3K36me3 = "#B0D9CB", EZH2 = "#6F8BA0",
  H3K27me3 = "#3E5B59", H3K9me3 = "#4A5F71")

## broad -> narrow (progressive zoom); base-resolution motif logo goes last
ZOOM_GROUPS <- list(
  "repressive / broad"   = list(half = 3000, label = "Repressive histone",
      marks = c("H3K36me3","EZH2","H3K27me3","H3K9me3")),
  "active histone"       = list(half = 1500, label = "Activating histone",
      marks = c("EP300","H3K27ac","H3K4me1","H3K4me3","RNA")),
  "CTCF / accessibility" = list(half = 500,  label = "CTCF / accessibility",
      marks = c("DNase","cohesin","CTCF","CAGE")))

## ENCODE SCREEN cCREs at the SPR locus (hg38), relative to the variant;
## variant lies inside EH38E2009446 (pELS, CTCF-bound).
ENCODE_CCRE <- tribble(
  ~name,          ~start, ~end,  ~type,
  "EH38E2009442",  -1432, -1169, "pELS",
  "EH38E2009443",  -1096,  -801, "pELS,CTCF-bound",
  "EH38E2009444",   -789,  -492, "PLS,CTCF-bound",
  "EH38E2009445",   -382,  -171, "pELS,CTCF-bound",
  "EH38E2009446",   -163,    70, "pELS,CTCF-bound",
  "EH38E2009447",    538,   867, "pELS,CTCF-bound",
  "EH38E2009448",    905,  1206, "pELS,CTCF-bound")
ENCODE_ELEMENT <- with(dplyr::filter(ENCODE_CCRE, start <= 0, end >= 0), c(start[1], end[1]))

## ---- data loading ----------------------------------------------------------
.spr_build_metadata <- function(bw_dir, mark_order = MARK_ORDER) {
  tibble(ref_path = list.files(bw_dir, pattern = "\\.ref\\.bw$",
                               recursive = TRUE, full.names = TRUE)) |>
    mutate(alt_path = str_replace(ref_path, "\\.ref\\.bw$", ".alt.bw"),
           mark     = str_remove(basename(ref_path), "\\.ref\\.bw$"),
           mark     = factor(mark, levels = mark_order)) |>
    arrange(mark) |> select(mark, ref_path, alt_path)
}
.spr_read_bw <- function(path, label) {
  as.data.frame(import(path)) |>
    transmute(score, type = factor(label, levels = c("reference", "alternate")),
              start)
}
.spr_load_track <- function(mark, meta, variant_pos = SPR_VARIANT_POS) {
  row <- filter(meta, mark == !!mark)
  bind_rows(.spr_read_bw(row$ref_path[1], "reference"),
            .spr_read_bw(row$alt_path[1], "alternate")) |>
    mutate(mark = factor(as.character(mark), levels = levels(meta$mark)),
           relative_pos = start - variant_pos)
}
.spr_load_all <- function(meta, variant_pos = SPR_VARIANT_POS)
  map_dfr(as.character(meta$mark), .spr_load_track, meta = meta, variant_pos = variant_pos)

## ---- % change at the regulatory element (each mark's ref peak +/- half) -----
pct_at_element <- function(df, search = c(-750, 750), half = 96) {
  d <- filter(df, relative_pos >= search[1], relative_pos <= search[2])
  peaks <- d |> filter(type == "reference") |>
    group_by(mark) |> slice_max(score, n = 1, with_ties = FALSE) |>
    transmute(mark, peak_pos = relative_pos)
  d |> inner_join(peaks, by = "mark") |>
    filter(relative_pos >= peak_pos - half, relative_pos <= peak_pos + half) |>
    group_by(mark, type) |> summarise(sig = mean(score), .groups = "drop") |>
    pivot_wider(names_from = type, values_from = sig) |>
    mutate(pct = 100 * (`alternate` - reference) / reference) |> arrange(mark)
}
.spr_labeller <- function(df, search = c(-750, 750), half = 96) {
  pc <- pct_at_element(df, search, half)
  as_labeller(setNames(sprintf("%s\n(%+.0f%%)", pc$mark, pc$pct), as.character(pc$mark)))
}

## ---- grouped browser (broad->narrow), cool per-mark fill + red alt ----------
plot_browser_grouped <- function(df, groups = ZOOM_GROUPS, mark_colors = MARK_COLORS,
                                 alt_col = "#d1495b", fill_alpha = 0.25, ref_lwd = 0.7,
                                 min_pct = 10, drop = c("RNA"),
                                 element = NULL, hl_fill = "#F4A261", hl_alpha = 0.22,
                                 bottom = NULL, bottom_height = 2.5) {
  pc <- pct_at_element(df)
  keep_marks <- setdiff(as.character(pc$mark[abs(pc$pct) >= min_pct]), drop)
  groups <- purrr::keep(map(groups, \(g) { g$marks <- intersect(g$marks, keep_marks); g }),
                        \(g) length(g$marks) > 0)
  df <- filter(df, mark %in% keep_marks)
  if (is.null(element)) element <- ENCODE_ELEMENT

  grp_panel <- function(marks, half, x_title = FALSE) {
    d <- df |> filter(mark %in% marks) |>
         mutate(mark = factor(as.character(mark), levels = marks)) |>
         filter(relative_pos >= -half, relative_pos <= half)
    fmt <- \(v) ifelse(v < 3, sprintf("%.1f", v), sprintf("%.0f", v))
    rng <- d |> group_by(mark) |>
           summarise(lab = sprintf("[0-%s]", fmt(max(score))), .groups = "drop")
    ggplot(d, aes(relative_pos, score)) +
      annotate("rect", xmin = element[1], xmax = element[2], ymin = -Inf, ymax = Inf,
               fill = hl_fill, alpha = hl_alpha) +
      geom_area(data = \(x) filter(x, type == "reference"), aes(fill = mark), alpha = fill_alpha) +
      geom_line(data = \(x) filter(x, type == "reference"), aes(colour = mark), linewidth = ref_lwd) +
      geom_line(data = \(x) filter(x, type == "alternate"), colour = alt_col, linewidth = .6) +
      geom_vline(xintercept = 0, colour = "grey30", linewidth = .4) +
      geom_text(data = rng, aes(x = -half, y = Inf, label = lab), inherit.aes = FALSE,
                hjust = -0.08, vjust = 1.35, size = 2.3, colour = "grey45") +
      scale_fill_manual(values = mark_colors, guide = "none") +
      scale_colour_manual(values = mark_colors, guide = "none") +
      scale_x_continuous(breaks = c(-half, -half/2, 0, half/2, half)) +
      facet_grid(rows = vars(mark), scales = "free_y", switch = "y",
                 labeller = .spr_labeller(df |> filter(mark %in% marks))) +
      labs(x = if (x_title) "position relative to variant (bp)" else NULL, y = NULL) +
      theme_classic(base_size = 10) +
      theme(strip.placement = "outside", strip.background = element_blank(),
            strip.text.y.left = element_text(angle = 0, size = 9.5), panel.spacing = unit(2, "pt"),
            axis.text.y = element_blank(), axis.ticks.y = element_blank(),
            axis.line.y = element_blank(), legend.position = "none")
  }
  panels <- map(groups, \(g) grp_panel(g$marks, g$half))
  rel_h  <- map_dbl(groups, \(g) length(g$marks) + 1)
  if (!is.null(bottom)) { panels <- c(panels, list(bottom)); rel_h <- c(rel_h, bottom_height) }
  body <- plot_grid(plotlist = panels, ncol = 1, align = "v", axis = "lr", rel_heights = rel_h)

  ctr_top <- cumsum(rel_h)[seq_along(groups)] - rel_h[seq_along(groups)] / 2
  yf <- 1 - ctr_top / sum(rel_h)
  right <- ggdraw()
  for (i in seq_along(groups))
    right <- right + draw_label(map_chr(groups, "label")[i], x = 0.55, y = yf[i],
                                angle = 270, size = 8.5, fontface = "bold", colour = "grey30")
  plot_grid(body, right, ncol = 2, rel_widths = c(1, 0.05))
}

## ---- base-resolution attribution motif logo --------------------------------
motif_logo <- function(attr_file, window = c(-8, 12), target = "CTCF", hl = "#d1495b",
                       hl_extend_right = 0) {
  a <- read_tsv(attr_file, show_col_types = FALSE) |> filter(rel >= window[1], rel <= window[2])
  bases <- c("A","C","G","T")
  mat <- matrix(0, 4, nrow(a), dimnames = list(bases, as.character(a$rel)))
  for (i in seq_len(nrow(a))) mat[a$base[i], i] <- a$attr[i] / 1000
  var_cols <- which(a$deleted)
  p <- ggseqlogo(mat, method = "custom")
  if (length(var_cols))
    p <- p + annotate("rect", xmin = min(var_cols) - 0.5, xmax = max(var_cols) + 0.5 + hl_extend_right,
                      ymin = -Inf, ymax = Inf, fill = hl, alpha = 0.18)
  p +
    annotate("text", x = 0.5, y = Inf, label = sprintf("[%.0f, %.0f]", min(mat), max(mat)),
             hjust = 0, vjust = 1.35, size = 2.3, colour = "grey45") +
    scale_x_continuous(breaks = seq(1, nrow(a), 4), labels = a$rel[seq(1, nrow(a), 4)]) +
    labs(x = "position relative to variant (bp)", y = sprintf("%s\nattribution", target)) +
    theme_classic(base_size = 10) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
          axis.title.y = element_text(size = 8, colour = "grey20", angle = 0, vjust = 0.5))
}

## ---- top-level builder -----------------------------------------------------
make_spr_browser <- function(bw_dir = SPR_BW_DIR, attr_file = SPR_ATTR) {
  meta       <- .spr_build_metadata(bw_dir)
  all_tracks <- .spr_load_all(meta)
  ## ref allele GCGCC spans 5 bp but only 4 are flagged deleted; extend one base right
  logo       <- motif_logo(attr_file, target = "CTCF", hl_extend_right = 1)
  plot_browser_grouped(all_tracks, bottom = logo, bottom_height = 3.6)
}
