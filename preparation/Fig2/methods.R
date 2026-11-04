# Data methods retained from the current local Fig2/00_setup.R.
# The source SHA-256 is recorded in the preparation manifest.
# Source result tables are the boundary; no association mapping is run.
suppressPackageStartupMessages(library(tidyverse))

add_derived <- function(df) {
  df %>% mutate(
    rare       = gvs_max_af < 0.01,
    odds_ratio = as.numeric(gvs_odds_ratio),
    p_value    = as.numeric(gvs_p_value)
  )
}

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
