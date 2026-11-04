# Scientific preparation functions retained from the current local Fig4/00_setup.R.
# Starts with GWAS-colocalization and conditional-TWAS result tables.
# Conditional TWAS and upstream association models are not rerun.
suppressPackageStartupMessages({library(tidyverse); library(data.table)})

RATE_COLUMNS <- c(stringent = "pct_stringent_union", lenient = "pct_lenient_union")

BLOOD_TRAITS <- c("eosinophil", "hematocrit", "lymphocyte", "MCH", "MCHC",
                  "monocyte", "neutrophil", "platelet_count", "RDW")

read_any <- function(path, ...) {
  args <- list(showProgress = FALSE, ...)
  if (grepl("\\.gz$", path)) args$cmd <- paste("gzip -dc", shQuote(path))
  else                         args$input <- path
  as_tibble(do.call(data.table::fread, args))
}

require_columns <- function(tbl, cols, label) {
  missing <- setdiff(cols, names(tbl))
  if (length(missing) > 0)
    stop(label, " is missing required column(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
}

as_logical_flag <- function(x) {
  if (is.logical(x)) return(replace_na(x, FALSE))
  if (is.numeric(x)) return(replace_na(x != 0, FALSE))
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

build_plot_data <- function(gene, coloc_rate, min_coloc_genes = 1L,
                            rate_type = "stringent") {
  rate_type <- tolower(trimws(rate_type))
  if (length(rate_type) != 1 || is.na(rate_type) || !rate_type %in% names(RATE_COLUMNS))
    stop("rate_type must be one of: ", paste(names(RATE_COLUMNS), collapse = ", "),
         call. = FALSE)
  rate_column <- unname(RATE_COLUMNS[[rate_type]])
  require_columns(gene, c("gene", "trait", "trait_category", "any_coloc"), "gene table")
  require_columns(coloc_rate, c("trait", "layer", unname(RATE_COLUMNS)), "coloc_rate table")
  if (length(min_coloc_genes) != 1 || is.na(min_coloc_genes) ||
      min_coloc_genes < 1 || min_coloc_genes != as.integer(min_coloc_genes))
    stop("min_coloc_genes must be a positive integer.", call. = FALSE)

  gene$any_coloc <- as_logical_flag(gene$any_coloc)
  coloc_counts <- gene %>%
    filter(any_coloc, !is.na(gene), gene != "", !is.na(trait_category), trait_category != "") %>%
    distinct(gene, trait, trait_category) %>%
    count(trait, trait_category, name = "count")

  pct_coloc <- coloc_rate %>%
    filter(layer == "union") %>%
    transmute(trait, pct_stringent_union, pct_lenient_union,
              pct_coloc_union = .data[[rate_column]])
  dup <- pct_coloc %>% count(trait) %>% filter(n > 1) %>% pull(trait)
  if (length(dup) > 0)
    stop("coloc_rate has multiple union rows for trait(s): ",
         paste(dup, collapse = ", "), call. = FALSE)

  combined <- coloc_counts %>% left_join(pct_coloc, by = "trait")
  missing_rate <- combined %>% filter(is.na(pct_coloc_union)) %>% pull(trait)
  if (length(missing_rate) > 0)
    message("Dropped (no union ", rate_type, " coloc rate): ",
            paste(unique(missing_rate), collapse = ", "))
  combined <- combined %>%
    filter(!is.na(pct_coloc_union)) %>% mutate(rate_type = .env$rate_type)
  trimmed <- combined %>% filter(count < min_coloc_genes) %>% pull(trait)
  if (length(trimmed) > 0)
    message("Trimmed ", length(trimmed), " trait(s) with fewer than ",
            min_coloc_genes, " colocalizing genes.")
  combined %>% filter(count >= min_coloc_genes)
}

coloc_by_trait <- function(refresh = FALSE)
  build_or_cache("coloc_by_trait", refresh = refresh, {
    build_plot_data(read_any(input_path("coloc_gene")),
                    read_any(input_path("coloc_rate")),
                    min_coloc_genes = 1, rate_type = "lenient")
  })

read_conditional_twas <- function(refresh = FALSE)
  build_or_cache("conditional_twas", refresh = refresh, {
    read_any(input_path("conditional_twas"))
  })

HEME_IMMUNE_KW <- str_c(c(
  "an[a]?emia", "thrombocyt", "neutropen", "leuk[a]?emia", "lymphoma",
  "h[a]?emoglobin", "platelet", "h[a]?emolytic", "spherocytosis", "sideroblastic",
  "myelodysplas", "granulocyt", "erythrocyt", "reticulocyt", "coagulation",
  "thromb", "bleeding", "h[a]?emophilia", "neutrophil", "eosinophil", "basophil",
  "monocyt", "myeloid", "h[a]?ematopoie", "pancytopenia", "h[a]?emochromatosis",
  "polycythemia", "thalass", "macrocyt", "blood cell", "bone marrow", "hemostasis",
  "porphyria", "fanconi", "immunodefic", "immune", "autoimmun", "autoinflammat",
  "immunodysreg", "agammaglobulin", "lymphoproliferat", "complement",
  "immunoglobulin", "hypogammaglobulin", "granulomatous", "histiocyt",
  "hemophagocytic", "interferonopathy"), collapse = "|")

derive_gene_sets <- function(clingen_df, drug_df) {
  clingen_cleaned <- clingen_df %>%
    distinct(targetId, targetFromSourceId, studyId, confidence,
             allelicRequirements, diseaseFromSource) %>%
    filter(confidence == "Definitive" | confidence == "Strong") %>%
    unnest(allelicRequirements)
  heme_immune_genes <- clingen_cleaned %>%
    filter(str_detect(tolower(diseaseFromSource), HEME_IMMUNE_KW)) %>%
    pull(targetId) %>% unique()

  dg <- as.data.table(drug_df)
  dg[, sym := toupper(gene_name)]
  druggable <- unique(dg[sym != "", sym])
  approved  <- unique(dg[approved == TRUE & sym != "", sym])

  list(druggable = druggable, approved = approved,
       mendelian = unique(clingen_cleaned$targetId),
       heme_immune = heme_immune_genes, clingen_cleaned = clingen_cleaned)
}

gene_sets <- function(refresh = FALSE)
  build_or_cache("gene_sets", refresh = refresh, {
    clingen_dir <- input_path("clingen")
    pq <- list.files(clingen_dir, pattern = "\\.parquet$", full.names = TRUE)
    if (length(pq) == 0) stop("No .parquet files under ", clingen_dir)
    clingen_df <- arrow::read_parquet(pq[1]) %>%
      select(targetFromSourceId, allelicRequirements, confidence, score,
             studyId, targetId, diseaseFromSource)
    drug_df <- read_any(input_path("drug_interactions"))
    derive_gene_sets(clingen_df, drug_df)
  })

prepare_qq_match_source <- function(data) {
  data %>%
    filter(status == "ok") %>%
    distinct(gene, coloc_trait, molecular_trait_id, .keep_all = TRUE) %>%
    mutate(any_coloc = coalesce(any_coloc, FALSE),
           inflation = abs(z_rare_cond) / pmax(abs(rare_zscore), 0.5),
           artifact  = coalesce(inflation > 3 & abs(Vcor) > 0.3, FALSE)) %>%
    filter(!artifact) %>%
    group_by(gene, coloc_trait) %>%
    arrange(desc(coalesce(GLCP, -Inf)), p_rare_cond, .by_group = TRUE) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    mutate(coloc = as.integer(any_coloc), log_tpm = log1p(mean_TPM),
           log_shet = log10(shet), log_n_rare = log1p(n_rare), abs_vcor = abs(Vcor)) %>%
    filter(is.finite(log_tpm), is.finite(log_shet),
           is.finite(log_n_rare), is.finite(abs_vcor))
}

calculate_qq_data <- function(matched_data, traits, common_cap = 40, rare_cap = 12.5) {
  matched_data %>%
    filter(as.character(coloc_trait) %in% traits) %>%
    select(group, common = common_pvalue, rare = p_rare_cond) %>%
    pivot_longer(c(common, rare), names_to = "test", values_to = "p_value") %>%
    filter(is.finite(p_value), between(p_value, 0, 1)) %>%
    group_by(test, group) %>%
    arrange(p_value, .by_group = TRUE) %>%
    mutate(rank = row_number(), n = n(),
           expected = -log10(rank / (n + 1)),
           observed = -log10(pmax(p_value, .Machine$double.xmin)),
           y_cap = case_when(test == "common" ~ common_cap,
                             test == "rare" ~ rare_cap, .default = Inf),
           observed_display = pmin(observed, y_cap),
           capped = observed > y_cap) %>%
    ungroup()
}

match_colocalized_pairs <- function(data, seed = 20260713, ratio = 1,
                                    caliper = 0.2, replace = FALSE) {
  set.seed(seed)
  match_fit <- MatchIt::matchit(
    coloc ~ log_tpm + log_shet + log_n_rare + abs_vcor, data = data,
    method = "nearest", distance = "glm", exact = ~ coloc_trait,
    ratio = ratio, replace = replace, caliper = caliper, std.caliper = TRUE)
  matched_data <- MatchIt::match.data(match_fit, drop.unmatched = TRUE) %>%
    mutate(group = if_else(coloc == 1L, "Colocalized genes", "Matched non-colocalized genes"),
           group = factor(group, levels = c("Matched non-colocalized genes", "Colocalized genes")))
  list(fit = match_fit, data = matched_data)
}

qq_match_source <- function(refresh = FALSE)
  build_or_cache("qq_match_source", refresh = refresh, {
    match_colocalized_pairs(prepare_qq_match_source(read_conditional_twas()))$data
  })

gene_twas_data <- function(refresh = FALSE)
  build_or_cache("gene_twas_data", refresh = refresh, {
    gs  <- gene_sets()
    adj <- read_conditional_twas() %>%
      distinct(gene, coloc_trait, molecular_trait_id, .keep_all = TRUE) %>%
      mutate(inflation = abs(z_rare_cond) / pmax(abs(rare_zscore), 0.5),
             artifact_primary = coalesce(inflation > 3 & abs(Vcor) > 0.3, FALSE)) %>%
      filter(!artifact_primary, any_coloc == TRUE) %>%
      group_by(TWAS_path) %>%
      mutate(fdr_rare_cond   = p.adjust(p_rare_cond,   method = "BH"),
             fdr_common      = p.adjust(common_pvalue, method = "BH"),
             rare_zscore_fdr = p.adjust(rare_pvalue,   method = "BH"),
             fdr_full        = p.adjust(full_pvalue,   method = "BH")) %>%
      ungroup() %>%
      mutate(genuine_rare = p_rare_cond < 0.05 & rare_pvalue < 0.05,
             t_fdrcond    = fdr_rare_cond < 0.10 & rare_pvalue < 0.05,
             t_fdrboth    = fdr_rare_cond < 0.10 & rare_zscore_fdr < 0.10,
             gene = factor(gene), coloc_trait = factor(coloc_trait),
             druggable = toupper(as.character(Description)) %in% gs$druggable,
             approved  = toupper(as.character(Description)) %in% gs$approved)

    # Guard: the artifact rule must have removed all offending rows.
    stopifnot(nrow(filter(adj, abs(Vcor) > 0.3, inflation > 3)) == 0)

    adj %>%
      group_by(gene, Description) %>%
      summarise(max_rare    = max(abs(z_rare_cond),   na.rm = TRUE),
                max_common  = max(abs(common_zscore), na.rm = TRUE),
                TPM         = mean(mean_TPM, na.rm = TRUE),
                druggable   = first(druggable), approved = first(approved),
                ever_common = any(fdr_common < 0.1, na.rm = TRUE),
                ever_rare   = any(t_fdrcond),
                shet        = first(shet), .groups = "drop") %>%
      mutate(log_shet = log10(shet), log_tpm = log10(TPM + 1),
             MendelianGene  = as.character(gene) %in% gs$clingen_cleaned$targetId,
             HemeImmuneGene = as.character(gene) %in% gs$heme_immune) %>%
      filter(!is.na(log_shet))
  })

fit_enrich <- function(outcome, predictors, data)
  glm(reformulate(predictors, response = outcome), data = data, family = binomial()) %>%
    broom::tidy(conf.int = TRUE, exponentiate = TRUE)

ENRICH_OUTCOMES <- c(
  approved       = "Approved Drug Target",
  druggable      = "Druggable gene",
  MendelianGene  = "Mendelian Gene\n(Any)",
  HemeImmuneGene = "Mendelian Gene\n(Immune/Hematological)")

ENRICH_PREDICTORS <- list(
  continuous = c("max_rare", "max_common", "log_shet", "log2(TPM)"),
  binary     = c("ever_rare", "ever_common", "log_shet", "log2(TPM)"))
