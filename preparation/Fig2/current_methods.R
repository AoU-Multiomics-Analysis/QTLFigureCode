# Current display-data methods from the November 2 figure sources.
# The input unit is a saved QTL result row. No participant-level input is used.

fig2_frequency_counts <- function(assay_tables) {
  dat <- purrr::imap_dfr(assay_tables, ~mutate(.x, Assay = .y))
  stopifnot(all(is.na(dat$pip) | between(dat$pip, 0, 1)))
  kept <- dat %>% filter(pip > 0.9, !is.na(gvs_max_af), gvs_max_af > 0)
  stopifnot(all(is.finite(kept$gvs_max_af)), all(kept$gvs_max_af <= 0.5))
  bins <- c("<1%", "1–<5%", "5–<10%", "10–50%")
  kept <- kept %>% mutate(Assay = factor(Assay, names(assay_tables)),
                         frequency_bin = factor(case_when(gvs_max_af < 0.01 ~ bins[1],
                                                          gvs_max_af < 0.05 ~ bins[2],
                                                          gvs_max_af < 0.10 ~ bins[3], TRUE ~ bins[4]), bins))
  counts <- kept %>% count(Assay, frequency_bin, .drop = FALSE) %>%
    group_by(Assay) %>% mutate(total = sum(n), percent = 100 * n / total) %>% ungroup()
  qc <- dat %>% group_by(Assay) %>%
    summarise(input_rows = n(), pip_above_09 = sum(pip > 0.9, na.rm = TRUE),
              missing_af = sum(pip > 0.9 & is.na(gvs_max_af), na.rm = TRUE),
              nonpositive_af = sum(pip > 0.9 & gvs_max_af <= 0, na.rm = TRUE), .groups = "drop")
  list(counts = counts, qc = qc)
}

fig2_ancestry_data <- function(assay_tables) {
  columns <- c("variant", "pip", "gvs_max_subpop", "gvs_max_af", "odds_ratio", "p_value", "rare")
  skew_df <- map_dfr(assay_tables, ~select(add_derived(.x), all_of(columns))) %>%
    filter(pip > 0.9, !is.na(gvs_max_subpop)) %>% distinct(variant, .keep_all = TRUE) %>%
    mutate(padj = p.adjust(p_value, method = "fdr"), enriched = odds_ratio > 5 & padj < 0.05,
           pop = toupper(gvs_max_subpop), rare_lab = if_else(rare, "Rare", "Common"))
  skew_counts <- skew_df %>% filter(enriched) %>% distinct(variant, pop) %>% count(pop)
  or_min <- min(skew_df$odds_ratio[skew_df$odds_ratio > 0], na.rm = TRUE)
  pair_line_df <- skew_df %>% distinct(pop) %>%
    mutate(x = as.numeric(factor(pop)) - 0.32, xend = as.numeric(factor(pop)) + 0.32,
           y = or_min * 0.55, yend = or_min * 0.55)
  label_df <- tibble(rare_lab = c("Common", "Rare"), x = 1 + c(-0.2, 0.2), y = or_min * 0.32)
  list(skew_df = skew_df, skew_counts = skew_counts, or_min = or_min,
       pair_line_df = pair_line_df, label_df = label_df)
}

# The saved enrichment result carries the four contingency counts. Its OR and
# Wald interval use 0.5 in every cell. Its two-sided Fisher P uses unmodified
# integer counts. These formulas reproduce the local source statistics.
fig2_outlier_statistics <- function(counts) {
  counts %>% mutate(
    focal_prop = focal_outliers / (focal_outliers + focal_non_outliers),
    bg_prop = bg_outliers / (bg_outliers + bg_non_outliers),
    odds_ratio = ((focal_outliers + 0.5) * (bg_non_outliers + 0.5)) /
      ((focal_non_outliers + 0.5) * (bg_outliers + 0.5)),
    log_odds_ratio = log(odds_ratio),
    se_log_or = sqrt(1 / (focal_outliers + 0.5) + 1 / (focal_non_outliers + 0.5) +
                       1 / (bg_outliers + 0.5) + 1 / (bg_non_outliers + 0.5)),
    ci_low = exp(log_odds_ratio - 1.96 * se_log_or),
    ci_high = exp(log_odds_ratio + 1.96 * se_log_or),
    p_value = pmap_dbl(pick(focal_outliers, focal_non_outliers, bg_outliers, bg_non_outliers),
                      function(...) fisher.test(matrix(c(...), nrow = 2, byrow = TRUE))$p.value))
}

fig2_median_data <- function(data) {
  data %>% filter(gene_type == "protein_coding") %>%
    mutate(pct_bin = if_else(as.character(pct_bin) == "(-10,-10)", "(-10,10)", as.character(pct_bin)),
           gene_type = "Protein Coding", lower = as.numeric(str_extract(pct_bin, "-?[0-9]+\\.?[0-9]*")))
}

fig2_cumulative_display <- function(cumulative) {
  plot_df <- cumulative$cumulative_df %>%
    mutate(direction = recode(direction, down = "Expression\ndecrease", up = "Expression\nincrease")) %>%
    filter(expected_het_expr < 4) %>% mutate(n_variants = pmin(n_variants, 4),
                                           label = if_else(n_variants < 4, as.character(n_variants), "4+"),
                                           k_label = factor(label, c("1", "2", "3", "4+"))) %>% distinct()
  labels <- plot_df %>% mutate(hit = if_else(direction == "Expression\ndecrease", expected_het_expr <= 0.6, expected_het_expr >= 1.5)) %>%
    group_by(direction, k_label) %>% summarise(n = n(), n_hit = sum(hit, na.rm = TRUE), prop = n_hit / n, .groups = "drop")
  list(plot_df = plot_df, label_df = labels)
}

# Annotation priority and regression are retained from Fig2/01_prep_data.R.
fig2_afc_data <- function(expression) {
  eQTL_COMB_susie <- add_derived(expression)
  remove_list <- c(
    "stop_lost", "inframe_deletion", "stop_retained_variant", "inframe_insertion",
    "CTCF.bound", "mature_mirna_variant", "upstream_gene_variant", "intron_variant",
    "intergenic_variant", "FANTOM5", "CTCF.only", "downstream_gene_variant",
    "synonymous_variant", "start_lost", "missense_variant",
    "non_coding_transcript_variant", "non_coding_transcript_exon_variant",
    "nmd_transcript_variant", "coding_sequence_variant", "NMD_transcript_variant")
  consequence_order <- c(
    "stop_gained", "frameshift_variant", "stop_lost", "start_lost",
    "splice_acceptor_variant", "splice_donor_variant", "splice_donor_5th_base_variant",
    "splice_donor_region_variant", "splice_polypyrimidine_tract_variant",
    "splice_region_variant", "PLS", "pELS", "dELS", "inframe_insertion",
    "inframe_deletion", "missense_variant", "coding_sequence_variant",
    "5_prime_UTR_variant", "3_prime_UTR_variant", "intron_variant",
    "upstream_gene_variant", "downstream_gene_variant", "intergenic_variant",
    "CTCF.bound", "CTCF.only", "DNase.H3K4me3", "FANTOM5")

  worst_csq <- eQTL_COMB_susie %>%
    filter(pip > 0.9, consequence != ".") %>%
    select(variant, log2_aFC, pip, rare, posterior_mean, consequence, molecular_trait_id,
           FANTOM5, pELS, dELS, PLS, `CTCF.bound`, `DNase.H3K4me3`, `CTCF.only`) %>%
    separate_rows(consequence, sep = ",") %>%
    mutate(value = TRUE) %>%
    pivot_wider(names_from = consequence, values_from = value, values_fill = FALSE) %>%
    pivot_longer(!c(variant, log2_aFC, pip, rare, posterior_mean, molecular_trait_id)) %>%
    mutate(name = trimws(name)) %>%
    filter(!name %in% remove_list, !str_detect(name, "DNase"), value == TRUE) %>%
    mutate(consequence_rank = match(name, consequence_order)) %>%
    filter(!is.na(consequence_rank)) %>%
    group_by(variant) %>% slice_min(consequence_rank, n = 1, with_ties = FALSE) %>% ungroup() %>%
    mutate(name = if_else(str_detect(name, "splice"), "Splice region", name))

  afc_plot_df <- worst_csq %>%
    mutate(abs_log2_aFC = pmin(abs(log2_aFC), 3),
           direction = if_else(log2_aFC >= 0, "Up", "Down")) %>%
    mutate(name = str_replace(name, "_prime_", "' "),
           name = str_replace(name, "UTR_variant", " UTR"),
           name = case_when(name == "stop_gained" ~ "Stop Gain",
                            name == "frameshift_variant" ~ "Frameshift", TRUE ~ name)) %>%
    mutate(name = fct_reorder(name, abs_log2_aFC, .fun = median, .desc = TRUE)) %>%
    mutate(ExpectedHet = (1 + 2^log2_aFC) / 2,
           high_effect = ExpectedHet < 0.6 | ExpectedHet > 1.4,
           color_group = paste0(name, "_", rare))

  log2_afc_regression <- afc_plot_df %>%
    group_by(name) %>% do(broom::tidy(lm(abs(log2_aFC) ~ rare, data = .))) %>%
    filter(!str_detect(term, "Int")) %>% ungroup() %>%
    mutate(padj = p.adjust(p.value, method = "fdr"),
           stars = case_when(padj < 0.001 ~ "***", padj < 0.01 ~ "**", padj < 0.05 ~ "*", TRUE ~ ""),
           label = paste0("+", round(estimate, 2), stars), x = 3.1)

  base_names  <- levels(afc_plot_df$name)
  base_cols   <- as.character(paletteer::paletteer_d("futurevisions::pegasi")[seq_along(base_names)])
  names(base_cols) <- base_names
  afc_palette <- c(setNames(colorspace::darken(base_cols, 0.3),  paste0(base_names, "_TRUE")),
                   setNames(colorspace::lighten(base_cols, 0.3), paste0(base_names, "_FALSE")))
  list(afc_plot_df = afc_plot_df, log2_afc_regression = log2_afc_regression, afc_palette = afc_palette)
}

# S16 keeps the original row-level BH universe; it does not deduplicate variants.
fig2_ancestry_bins <- function(expression) {
  eQTL_COMB_susie <- add_derived(expression)
  af_cols <- c("gvs_all_af", "gvs_max_af", "gvs_afr_af", "gvs_amr_af",
               "gvs_eas_af", "gvs_eur_af", "gvs_mid_af", "gvs_sas_af", "gvs_oth_af")

  AncestrySkewAnnotated <- eQTL_COMB_susie %>%
    filter(pip > 0.9) %>%                        # fine-mapped variants only
    mutate(across(all_of(af_cols), as.numeric)) %>%
    mutate(
      MAF     = if_else(gvs_all_af > 0.5, 1 - gvs_all_af, gvs_all_af),
      log_MAF = -log10(MAF),
      padj    = p.adjust(p_value, method = "fdr"),
      af_bin  = ntile(log_MAF, 30),
      enriched = odds_ratio > 5 & padj < 0.05,
      present_afr = gvs_afr_af > 0, present_eur = gvs_eur_af > 0,
      present_eas = gvs_eas_af > 0, present_amr = gvs_amr_af > 0,
      present_sas = gvs_sas_af > 0, present_mid = gvs_mid_af > 0,
      present_oth = gvs_oth_af > 0
    ) %>%
    mutate(shared_all_continents = present_afr & present_eur & present_eas &
             present_amr & present_mid & present_sas & present_oth)

  # AF bins -> representative x (geometric mean of the bin's AF range).
  bin_key <- AncestrySkewAnnotated %>%
    group_by(af_bin) %>%
    summarise(min_log = min(log_MAF, na.rm = TRUE),
              max_log = max(log_MAF, na.rm = TRUE), .groups = "drop") %>%
    mutate(min_af = 10^(-max_log), max_af = 10^(-min_log), x = sqrt(min_af * max_af)) %>%
    select(af_bin, x)

  all_pops   <- sort(unique(AncestrySkewAnnotated$gvs_max_subpop))
  bin_totals <- AncestrySkewAnnotated %>% count(af_bin, name = "total_variants")

  # Enriched-variant composition by ancestry, per AF bin.
  stack_df <- AncestrySkewAnnotated %>%
    filter(enriched) %>%
    count(af_bin, gvs_max_subpop, name = "n_enriched") %>%
    complete(af_bin = bin_key$af_bin, gvs_max_subpop = all_pops, fill = list(n_enriched = 0)) %>%
    left_join(bin_totals, by = "af_bin") %>%
    mutate(total_variants = replace_na(total_variants, 0),
           prop = if_else(total_variants > 0, n_enriched / total_variants, 0)) %>%
    left_join(bin_key, by = "af_bin")

  # Proportion of variants shared across all continents, per AF bin (smoothed).
  line_df <- AncestrySkewAnnotated %>%
    group_by(af_bin) %>%
    summarise(prop_shared = mean(shared_all_continents, na.rm = TRUE), .groups = "drop") %>%
    complete(af_bin = bin_key$af_bin, fill = list(prop_shared = 0)) %>%
    left_join(bin_key, by = "af_bin") %>%
    arrange(x) %>%
    mutate(prop_shared_smooth = slider::slide_dbl(prop_shared, mean, .before = 2, .after = 2, .complete = FALSE))
  list(stack_df = stack_df, line_df = line_df)
}
