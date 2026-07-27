# =============================================================================
# Figure 4 — shared setup (sourced by every panel / assembly script).
# Run all scripts from the PROJECT ROOT. Each script starts with:
#   source("fig4/00_setup.R")
#
# Figure 4 covers colocalization of molecular QTLs with GWAS traits and the
# downstream rare-variant / conditional-TWAS analyses.
#
# See docs/superpowers/specs/2026-07-24-fig4-cleanup-design.md
# =============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(patchwork)
  library(ggh4x)
  library(MatchIt)
  library(arrow)
  library(broom)
  library(scales)
  library(showtext)
  library(sysfonts)
  library(paletteer)
})

# --- Paths (override via env vars if needed) --------------------------------
INPUT    <- Sys.getenv("FIG4_INPUT",    ".")            # staging root (cwd)
DATA_DIR <- Sys.getenv("FIG4_DATA",     "fig4/data")    # memoized derived tables
OBJ_DIR  <- Sys.getenv("FIG4_OBJ",      "fig4/objects") # per-panel ggplot objects
OUTDIR   <- Sys.getenv("FIG4_OUTDIR",   ".")            # final Figure4.pdf/png
MANIFEST <- Sys.getenv("FIG4_MANIFEST", "fig4/inputs.tsv")
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OBJ_DIR,  showWarnings = FALSE, recursive = TRUE)

source("common/manifest.R")

# --- Font: Inter via Google if reachable, else default sans -----------------
FONT <- "sans"
if (!"inter" %in% sysfonts::font_families()) {
  if (isTRUE(tryCatch({ font_add_google("Inter", "inter"); TRUE }, error = function(e) FALSE)))
    FONT <- "inter"
} else FONT <- "inter"

DPI <- 400                     # showtext DPI; must match every ggsave() DPI
showtext_opts(dpi = DPI)
showtext_auto()

# --- Design system ----------------------------------------------------------
# element_text() sizes are in pt; geom_text()/geom_label()/annotate() are in mm.
SZ_AXIS_TITLE <- 14; SZ_AXIS_TEXT <- 12; SZ_STRIP <- 12
SZ_LEG_TEXT   <- 10; SZ_LEG_TITLE <- 11; SZ_TAG   <- 14
SZ_TITLE      <- 12

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

# --- Trait-category palette / labels (used by Panel A) ----------------------
GROUP_LABELS <- c(
  autoimmune = "AutoImm", inflammatory_biomarker = "Inflam",
  allergic_respiratory = "Allerg", dermatological = "Derm",
  cardiometabolic = "CardioM", lipid_metabolic = "Lipid", lipid = "Lipid",
  metabolic = "Metabolic", clinical_chemistry = "ClinChem",
  hematological = "Heme", anthropometric = "Anthro", renal = "Renal",
  renal_urological = "Renal/Urol", gastrointestinal = "GI",
  reproductive_endocrine = "ReproEnd", musculoskeletal = "MSK",
  cancer = "Cancer", neurodegenerative = "Neurodeg", neurological = "Neurol",
  neuropsychiatric = "NeuroPsy", ophthalmological = "Ophth"
)
GROUP_COLORS <- c(
  autoimmune = "#D62828", inflammatory_biomarker = "#F77F00",
  allergic_respiratory = "#2A9D8F", dermatological = "#E76F51",
  cardiometabolic = "#277DA1", lipid_metabolic = "#6A994E", lipid = "#6A994E",
  metabolic = "#80B918", clinical_chemistry = "#F4A261",
  hematological = "#A23B72", anthropometric = "#E9C46A", renal = "#457B9D",
  renal_urological = "#8D99AE", gastrointestinal = "#BC6C25",
  reproductive_endocrine = "#E5989B", musculoskeletal = "#7F5539",
  cancer = "#9B5DE5", neurodegenerative = "#4361EE", neurological = "#3A0CA3",
  neuropsychiatric = "#B5179E", ophthalmological = "#48CAE4"
)
RATE_COLUMNS <- c(stringent = "pct_stringent_union", lenient = "pct_lenient_union")

# --- Blood traits used by Panel C -------------------------------------------
BLOOD_TRAITS <- c("eosinophil", "hematocrit", "lymphocyte", "MCH", "MCHC",
                  "monocyte", "neutrophil", "platelet_count", "RDW")

# =============================================================================
# Derived-data cache (identical contract to Fig3): lazy `expr`, cache to
# DATA_DIR/<name>.rds, keyed by NAME only. Clear DATA_DIR after changing an input.
# =============================================================================
build_or_cache <- function(name, expr, refresh = FALSE) {
  f <- file.path(DATA_DIR, paste0(name, ".rds"))
  if (!refresh && file.exists(f)) return(readRDS(f))
  message("building ", name, " ...")
  val <- force(expr)
  saveRDS(val, f)
  val
}

# =============================================================================
# Panel A helpers — coloc-by-trait-category plot
# =============================================================================
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

complete_category_colors <- function(categories) {
  categories <- sort(unique(as.character(categories)))
  missing <- setdiff(categories, names(GROUP_COLORS))
  fallback <- if (length(missing) > 0)
    setNames(scales::hue_pal(l = 55, c = 90)(length(missing)), missing) else character()
  c(GROUP_COLORS, fallback)
}

blend_with_white <- function(hex, amount = 0.65) {
  x <- grDevices::col2rgb(hex) / 255
  y <- x * (1 - amount) + amount
  grDevices::rgb(y[1, ], y[2, ], y[3, ])
}

make_fill <- function(category, intensity, colors) {
  base  <- unname(colors[as.character(category)])
  light <- blend_with_white(base, 0.72)
  vapply(seq_along(base), function(i) {
    palette <- grDevices::colorRampPalette(c(light[[i]], base[[i]]))(100)
    palette[pmax(1, pmin(100, round(intensity[[i]] * 99) + 1))]
  }, character(1))
}

reorder_within <- function(x, by, within, fun = mean, sep = "___")
  stats::reorder(paste(x, within, sep = sep), by, FUN = fun)

scale_y_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  scale_y_discrete(..., labels = function(x) gsub(reg, "", x))
}

# Build the per-trait count + union-coloc-rate table. `rate_type` selects which
# union rate column drives pct_coloc_union. Error messages name the argument
# plainly (the original used optparse-style "--gene" labels).
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

build_empty_plot <- function() {
  ggplot() +
    annotate("text", x = 0, y = 0,
             label = "No traits have both colocalizing genes and the selected union coloc rate",
             color = "grey35", size = 4) +
    xlim(-1, 1) + ylim(-1, 1) + theme_void()
}

build_plot <- function(combined) {
  colors <- complete_category_colors(combined$trait_category)
  plotdat <- combined %>%
    group_by(trait_category) %>%
    mutate(
      intensity = if (max(count) == min(count)) rep(0.75, n())
                  else (count - min(count)) / (max(count) - min(count)),
      intensity = 0.25 + 0.75 * intensity,
      fill_col  = make_fill(trait_category, intensity, colors)) %>%
    ungroup() %>%
    mutate(
      trait_category = forcats::fct_drop(forcats::fct_reorder(trait_category, count, .fun = sum)),
      trait          = reorder_within(trait, count, trait_category),
      base_col       = unname(colors[as.character(trait_category)]))

  long <- bind_rows(
    plotdat %>% transmute(trait_category, trait, base_col, fill_col, metric = "a_count", value = count),
    plotdat %>% transmute(trait_category, trait, base_col, fill_col, metric = "b_rate",  value = pct_coloc_union))
  xline <- -max(plotdat$count) * 0.03
  seg  <- plotdat %>% count(trait_category, name = "nbar") %>% mutate(metric = "a_count")
  vref <- tibble(metric = "b_rate", xi = median(plotdat$pct_coloc_union, na.rm = TRUE))
  metric_labels <- c(a_count = "Colocalizing Genes", b_rate = "CS coloc rate ")
  category_labeller <- c(GROUP_LABELS, setNames(
    setdiff(as.character(unique(plotdat$trait_category)), names(GROUP_LABELS)),
    setdiff(as.character(unique(plotdat$trait_category)), names(GROUP_LABELS))))

  ggplot(long, aes(y = trait)) +
    geom_col(data = ~ subset(.x, metric == "a_count"),
             aes(x = value, fill = fill_col), width = 0.85, show.legend = FALSE) +
    geom_segment(data = seg, inherit.aes = FALSE,
                 aes(y = 0.5, yend = nbar + 0.5, x = xline, xend = xline),
                 color = "grey60", linewidth = 0.5) +
    geom_vline(data = vref, aes(xintercept = xi),
               linetype = "dashed", color = "grey70", linewidth = 0.4) +
    geom_segment(data = ~ subset(.x, metric == "b_rate"),
                 aes(x = 0, xend = value, yend = trait, color = base_col), linewidth = 0.4) +
    geom_point(data = ~ subset(.x, metric == "b_rate"),
               aes(x = value, color = base_col), size = 1.6) +
    scale_fill_identity() + scale_color_identity() + scale_y_reordered() +
    facet_grid(trait_category ~ metric, scales = "free", space = "free_y", switch = "y",
               labeller = labeller(
                 trait_category = as_labeller(category_labeller, default = label_value),
                 metric = as_labeller(metric_labels))) +
    ggh4x::facetted_pos_scales(x = list(
      scale_x_continuous(),
      scale_x_continuous(labels = function(x) paste0(x, "%"),
                         expand = expansion(mult = c(0, 0.08))))) +
    ggh4x::force_panelsizes(cols = c(3, 1.3)) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_classic() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          panel.spacing.x = grid::unit(8, "pt"), panel.spacing.y = grid::unit(2, "pt"),
          strip.placement = "outside", strip.background = element_blank(),
          strip.text.y.left = element_text(angle = 0, hjust = 1, size = 12),
          strip.text.x = element_text(size = 12))
}

# Builder: the per-trait coloc table (lenient union rate, matching the figure).
coloc_by_trait <- function(refresh = FALSE)
  build_or_cache("coloc_by_trait", refresh = refresh, {
    build_plot_data(read_any(input_path("coloc_gene")),
                    read_any(input_path("coloc_rate")),
                    min_coloc_genes = 1, rate_type = "lenient")
  })

# =============================================================================
# Panel C/D upstream — conditional TWAS read, and gene-set annotations
# =============================================================================
read_conditional_twas <- function(refresh = FALSE)
  build_or_cache("conditional_twas", refresh = refresh, {
    read_any(input_path("conditional_twas"))
  })

# Heme / immune disease keyword pattern (hand-curated; see NOTES.md).
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

# Pure gene-set derivation (testable without arrow / the parquet read).
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

# =============================================================================
# Panel C helpers — conditional-TWAS QQ with propensity matching
# =============================================================================
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

plot_qq <- function(qq_data) {
  ggplot(qq_data, aes(x = expected, y = observed_display, colour = group)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                linewidth = 0.55, colour = "#666666") +
    geom_point(size = 1, alpha = 0.8) +
    facet_wrap(~test, scales = "free_y",
               labeller = as_labeller(c(common = "Common", rare = "Rare"))) +
    scale_colour_manual(values = c("Matched non-colocalized genes" = "#2B2B2B",
                                   "Colocalized genes" = "#275D8E"), name = NULL) +
    guides(colour = guide_legend(override.aes = list(shape = 16, size = 2.5, alpha = 1))) +
    labs(x = expression(Expected ~ -log[10](P)), y = expression(Observed ~ -log[10](P))) +
    theme_classic(base_size = 12) +
    theme(strip.background = element_blank(), strip.text = element_text(size = 12),
          panel.spacing = grid::unit(2, "lines"), legend.position = "top",
          legend.text = element_text(size = 10), legend.key.width = grid::unit(1.2, "lines"))
}

qq_match_source <- function(refresh = FALSE)
  build_or_cache("qq_match_source", refresh = refresh, {
    match_colocalized_pairs(prepare_qq_match_source(read_conditional_twas()))$data
  })

# =============================================================================
# Panel D upstream — colocalized conditional-TWAS enrichment
# -----------------------------------------------------------------------------
# NOTE: the original computed ColocAdjustedConditionalTWAS TWICE. The first copy
# referenced an undefined `fdr_rare_sig` and fed GeneTWASData; the second copy
# recomputed it cleanly (BH adjustment) but went unused. We keep ONLY the clean
# version and build GeneTWASData from it (see NOTES.md).
# =============================================================================
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
