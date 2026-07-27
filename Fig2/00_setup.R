# =============================================================================
# Figure 2 — shared setup (sourced by every prep / panel / assembly script)
# -----------------------------------------------------------------------------
# Run all scripts from the PROJECT ROOT (the folder containing input_data/ and
# fig2/). Each script starts with:  source("fig2/00_setup.R")
# Re-sourcing is cheap and idempotent.
# =============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(scales)
  library(patchwork)
  library(ggbeeswarm)
  library(showtext)
  library(sysfonts)
  library(paletteer)
  library(ggdist)
  library(ggh4x)
  library(rlang)
  library(slider)
  library(rcartocolor)
  library(colorspace)
})

# --- Paths (override via env vars if needed) --------------------------------
INPUT    <- Sys.getenv("FIG2_INPUT",    ".")            # staging root (cwd); inputs live under input_data/
DATA_DIR <- Sys.getenv("FIG2_DATA",     "fig2/data")    # per-panel data (01_prep)
OBJ_DIR  <- Sys.getenv("FIG2_OBJ",      "fig2/objects") # per-panel ggplot objects
OUTDIR   <- Sys.getenv("FIG2_OUTDIR", "fig2")         # final Fig2.pdf/png
MANIFEST <- Sys.getenv("FIG2_MANIFEST", "fig2/inputs.tsv")
for (d in c(DATA_DIR, OBJ_DIR, OUTDIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

source("common/manifest.R")

# --- Font: Inter via Google if reachable, else default sans -----------------
FONT <- "sans"
if (!"inter" %in% sysfonts::font_families()) {
  if (isTRUE(tryCatch({ font_add_google("Inter", "inter"); TRUE }, error = function(e) FALSE)))
    FONT <- "inter"
} else FONT <- "inter"

DPI <- 300                     # showtext DPI; must match every ggsave() DPI
showtext_opts(dpi = DPI)
showtext_auto()

# --- Design system: uniform text sizes --------------------------------------
# element_text() sizes are in pt; geom_text()/annotate() sizes are in mm.
SZ_AXIS_TITLE <- 22; SZ_AXIS_TEXT <- 18; SZ_STRIP <- 18
SZ_LEG_TEXT   <- 15; SZ_LEG_TITLE <- 17; SZ_TAG   <- 26
ANN_LG <- 8.0; ANN_MD <- 6.0; ANN_SM <- 4.5

theme_set(theme_classic())
theme_update(
  text             = element_text(family = FONT),
  axis.title       = element_text(size = SZ_AXIS_TITLE),
  axis.text        = element_text(size = SZ_AXIS_TEXT),
  strip.text       = element_text(size = SZ_STRIP),
  strip.background = element_blank(),
  legend.text      = element_text(size = SZ_LEG_TEXT),
  legend.title     = element_text(size = SZ_LEG_TITLE)
)

blank_plot <- ggplot() + theme_void()

# --- Shared helpers ---------------------------------------------------------
# Prettify GVS sub-population codes (afr -> AFR, oth -> ADMIX).
pretty_pop <- function(x) toupper(if_else(x == "oth", "admix", x))

# Derived columns the merged files no longer ship with.
add_derived <- function(df) {
  df %>% mutate(
    rare       = gvs_max_af < 0.01,
    odds_ratio = as.numeric(gvs_odds_ratio),
    p_value    = as.numeric(gvs_p_value)
  )
}

# Read a manifest key's file (single-file keys) with fread; the manifest maps
# the key -> input_data/<file> and stages it from gs_uri if ever added.
read_qtl <- function(key) fread(input_path(key), nThread = 4)

# For each gene, accumulate up/down variants by descending |effect| and track the
# implied heterozygous expression, then report how many genes reach dosage-like
# thresholds after including k variants. (Used by Panel D.)
make_cumulative_afc_curve_both <- function(df,
                                           gene_col = "gene",
                                           afc_col  = "log2_aFC",
                                           down_thresholds = c(0.75, 0.65, 0.60, 0.55),
                                           up_thresholds   = c(1.10, 1.25, 1.40, 1.50, 2.00, 2.50, 3.00),
                                           max_k = NULL) {
  gene_sym <- rlang::sym(gene_col)
  afc_sym  <- rlang::sym(afc_col)

  df_clean  <- df %>% dplyr::filter(!is.na(!!gene_sym), !is.na(!!afc_sym))
  all_genes <- df_clean %>% dplyr::distinct(!!gene_sym) %>% dplyr::transmute(gene = !!gene_sym)

  build_direction <- function(data, direction, thresholds) {
    df_use <- if (direction == "down") {
      data %>% dplyr::filter(!!afc_sym < 0) %>% dplyr::arrange(!!gene_sym, !!afc_sym)
    } else {
      data %>% dplyr::filter(!!afc_sym > 0) %>% dplyr::arrange(!!gene_sym, dplyr::desc(!!afc_sym))
    }

    if (nrow(df_use) == 0) {
      cumulative_df <- tibble::tibble(gene = character(), n_variants = integer(),
        cumulative_log2_aFC = numeric(), cumulative_aFC = numeric(),
        expected_het_expr = numeric(), direction = character())
      min_n_by_thresh <- tidyr::expand_grid(gene = all_genes$gene, threshold = thresholds) %>%
        dplyr::mutate(min_n = Inf, direction = direction)
      return(list(cumulative_df = cumulative_df, min_n_by_thresh = min_n_by_thresh))
    }

    cumulative_df <- df_use %>%
      dplyr::group_by(!!gene_sym) %>%
      dplyr::mutate(
        n_variants          = dplyr::row_number(),
        cumulative_log2_aFC = cumsum(!!afc_sym),
        cumulative_aFC      = 2^cumulative_log2_aFC,
        expected_het_expr   = (1 + cumulative_aFC) / 2
      ) %>%
      dplyr::ungroup() %>%
      dplyr::transmute(gene = !!gene_sym, n_variants, cumulative_log2_aFC,
                       cumulative_aFC, expected_het_expr, direction = direction)

    min_n_by_thresh <- purrr::map_dfr(thresholds, function(thresh) {
      hits <- if (direction == "down") {
        cumulative_df %>% dplyr::filter(expected_het_expr <= thresh)
      } else {
        cumulative_df %>% dplyr::filter(expected_het_expr >= thresh)
      }
      min_n <- hits %>% dplyr::group_by(gene) %>%
        dplyr::summarise(min_n = min(n_variants), .groups = "drop")
      all_genes %>%
        dplyr::left_join(min_n, by = "gene") %>%
        dplyr::mutate(threshold = thresh,
                      min_n = dplyr::if_else(is.na(min_n), Inf, as.numeric(min_n)),
                      direction = direction)
    })

    list(cumulative_df = cumulative_df, min_n_by_thresh = min_n_by_thresh)
  }

  down_res <- build_direction(df_clean, "down", down_thresholds)
  up_res   <- build_direction(df_clean, "up",   up_thresholds)

  cumulative_df   <- dplyr::bind_rows(down_res$cumulative_df, up_res$cumulative_df)
  min_n_by_thresh <- dplyr::bind_rows(down_res$min_n_by_thresh, up_res$min_n_by_thresh)

  if (is.null(max_k)) {
    finite_vals <- min_n_by_thresh$min_n[is.finite(min_n_by_thresh$min_n)]
    max_k <- if (length(finite_vals) == 0) 1 else max(finite_vals)
  }

  curve_df <- tidyr::expand_grid(direction = c("down", "up"),
                                 threshold = c(down_thresholds, up_thresholds),
                                 k = seq_len(max_k)) %>%
    dplyr::filter((direction == "down" & threshold %in% down_thresholds) |
                    (direction == "up" & threshold %in% up_thresholds)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(prop_genes = mean(
      min_n_by_thresh$min_n[min_n_by_thresh$direction == direction &
                              min_n_by_thresh$threshold == threshold] <= k)) %>%
    dplyr::ungroup()

  list(cumulative_df = cumulative_df, min_n_by_thresh = min_n_by_thresh, curve_df = curve_df)
}

# Panel D plot builder (kept here so both the panel script and any standalone
# use share one definition).
PlotExpectedExpression <- function(cumulative_effects_df) {
  decrease_band_high <- 0.6; decrease_band_low <- 0.45
  increase_band_high <- 4.0; increase_band_low <- 1.45

  band_df <- tibble::tibble(
    direction = c("Expression\ndecrease", "Expression\nincrease"),
    xmin = -Inf, xmax = Inf,
    ymin = c(decrease_band_low, increase_band_low),
    ymax = c(decrease_band_high, increase_band_high))
  band_label_df <- tibble::tibble(
    direction = c("Expression\ndecrease", "Expression\nincrease"),
    x = c(0.5, 4.05),
    y = c(decrease_band_low - 0.02, increase_band_high + 0.15),
    label = c("Deletion-like", "Duplication-like"), hjust = c(0, 0))

  plot_df <- cumulative_effects_df$cumulative_df %>%
    mutate(direction = recode(direction, down = "Expression\ndecrease",
                                          up   = "Expression\nincrease")) %>%
    filter(expected_het_expr < 4) %>%
    mutate(n_variants = pmin(n_variants, 4),
           label   = if_else(n_variants < 4, as.character(n_variants), "4+"),
           k_label = factor(label, levels = c("1", "2", "3", "4+"))) %>%
    distinct()

  label_df <- plot_df %>%
    mutate(hit = case_when(
      direction == "Expression\ndecrease" ~ expected_het_expr <= 0.6,
      direction == "Expression\nincrease" ~ expected_het_expr >= 1.5)) %>%
    group_by(direction, k_label) %>%
    summarise(n = dplyr::n(), n_hit = sum(hit, na.rm = TRUE),
              prop = n_hit / n, label = scales::percent(prop, accuracy = 1),
              .groups = "drop") %>%
    mutate(y = if_else(direction == "Expression\ndecrease", 1.1, 4.65))

  p <- ggplot(plot_df, aes(factor(k_label), expected_het_expr,
                           color = direction, fill = direction)) +
    geom_violin(width = 0.8, alpha = 0.3, color = NA, trim = FALSE, show.legend = FALSE) +
    geom_quasirandom(width = 0.1, size = 0.2, alpha = 0.2, show.legend = FALSE) +
    geom_hline(yintercept = 1, linetype = 2, color = "grey40") +
    facet_wrap(~direction, scales = "free_y", nrow = 2, strip.position = "top") +
    scale_fill_manual(values  = c("#B65C73", "#2A9D8F")) +
    scale_color_manual(values = c("#B65C73", "#2A9D8F")) +
    stat_summary(fun = median, geom = "crossbar", linewidth = 0.3, width = 0.2,
                 alpha = 0.2, color = "black", show.legend = FALSE) +
    ggh4x::facetted_pos_scales(y = list(
      scale_y_continuous(breaks = c(0.5, 0.7, 0.9), labels = percent_format(accuracy = 1)),
      scale_y_continuous(breaks = c(1.5, 2.5, 3.5), labels = percent_format(accuracy = 1)))) +
    geom_text(data = label_df, aes(k_label, y, label = label),
              inherit.aes = FALSE, size = ANN_LG) +
    geom_rect(data = band_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE, fill = "grey40", alpha = 0.05) +
    geom_text(data = band_label_df, aes(x, y, label = label, hjust = hjust),
              inherit.aes = FALSE, size = ANN_SM, color = "grey35", fontface = "bold") +
    labs(x = "Number of independent variants included",
         y = "Expected heterozygous expression (%)") +
    coord_cartesian(clip = "off") +
    theme(legend.position = "none",
          axis.line.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title.y.left = element_text(margin = margin(r = 30)),
          plot.margin = margin(r = 0, t = 50, b = 1, l = 100))

  p$facet$params$labeller <- as_labeller(c(
    "Expression\ndecrease" = "Number Decreasing Variants",
    "Expression\nincrease" = "Number Increasing Variants"))
  p
}
