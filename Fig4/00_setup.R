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
