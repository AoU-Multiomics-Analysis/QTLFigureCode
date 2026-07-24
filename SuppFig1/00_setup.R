# =============================================================================
# Supplemental Figure 1 — shared setup (sourced by every panel / extra / assembly
# script). Run all scripts from the PROJECT ROOT:
#   source("suppfig1/00_setup.R")
#
# SuppFig1 compares All of Us bulk RNA-seq expression to GTEx whole blood across
# gene types (coding correlation + detection, GO enrichment, TF expression,
# lncRNA), with additional per-gene-type comparisons in suppfig1/extra/.
#
# See docs/superpowers/specs/2026-07-24-suppfig1-cleanup-design.md
# =============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(rtracklayer)
  library(enrichR)
  library(hexbin)
  library(ggrepel)
  library(rcartocolor)
  library(broom)
  library(patchwork)
  library(showtext)
  library(sysfonts)
  library(paletteer)
})

# --- Paths (override via env vars if needed) --------------------------------
INPUT    <- Sys.getenv("SUPPFIG1_INPUT",    ".")                 # staging root (cwd)
DATA_DIR <- Sys.getenv("SUPPFIG1_DATA",     "suppfig1/data")     # memoized tables
OBJ_DIR  <- Sys.getenv("SUPPFIG1_OBJ",      "suppfig1/objects")  # per-panel objects
OUTDIR   <- Sys.getenv("SUPPFIG1_OUTDIR",   ".")                 # final outputs
MANIFEST <- Sys.getenv("SUPPFIG1_MANIFEST", "suppfig1/inputs.tsv")
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OBJ_DIR,  showWarnings = FALSE, recursive = TRUE)

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

# --- Design system ----------------------------------------------------------
# element_text() sizes are in pt; geom_text()/annotate() are in mm. This figure
# uses bordered panels (folded in from the original's theme_update).
SZ_AXIS_TITLE <- 20; SZ_AXIS_TEXT <- 15; SZ_STRIP <- 20
SZ_LEG_TEXT   <- 12; SZ_LEG_TITLE <- 14; SZ_TAG   <- 25

theme_set(theme_classic())
theme_update(
  text             = element_text(family = FONT),
  panel.border     = element_rect(fill = NA, linewidth = 0.9),
  axis.line        = element_blank(),
  axis.title       = element_text(size = SZ_AXIS_TITLE),
  axis.text        = element_text(size = SZ_AXIS_TEXT),
  strip.text       = element_text(size = SZ_STRIP),
  strip.background = element_blank(),
  legend.text      = element_text(size = SZ_LEG_TEXT),
  legend.title     = element_text(size = SZ_LEG_TITLE)
)

# =============================================================================
# Derived-data cache (identical contract to Fig3/Fig4): lazy `expr`, cache to
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
# Folded helpers
# =============================================================================
# Mean TPM per gene from a GCT-style matrix (Name, Description, then samples),
# joined to gene-type annotation. The original also emitted SE = sd/n, which is
# neither a standard error (sd/sqrt(n)) nor plotted anywhere; it is dropped.
compute_mean_TPM <- function(TPM_df, gene_annotations) {
  annotation_dat <- gene_annotations %>% distinct()
  TPM_df %>%
    pivot_longer(!c(Name, Description)) %>%
    group_by(Description, Name) %>%
    summarize(mean_TPM = mean(value), .groups = "drop") %>%
    left_join(annotation_dat, by = c("Name" = "gene_id"))
}

# Correlation scatter (AoU vs GTEx log2 TPM) beside cohort-specific detection
# counts. `hex = TRUE` uses geom_hex; after_stat() replaces the deprecated stat().
plot_gene_comparison <- function(TPM_df, point_size = 3, detection_threshold = 1,
                                 gene_type = "", hex = FALSE) {
  axis_label <- paste0("Unique ", gene_type, " genes")
  cor_values <- cor(log2(TPM_df$mean_TPM_AoU + 1), log2(TPM_df$mean_TPM_GTEx + 1),
                    use = "complete.obs")
  plot_text <- paste0("Cor = ", round(cor_values, digits = 2))

  cor_plot <- TPM_df %>%
    ggplot(aes(x = log2(mean_TPM_AoU + 1), y = log2(mean_TPM_GTEx + 1))) +
    xlab("log2(AoU mean TPM)") + ylab("log2(GTEx mean TPM)") +
    annotate("text", x = 1.8, y = 19, label = plot_text, size = 7)
  if (!hex) {
    cor_plot <- cor_plot +
      # point_size is honored here; the original hardcoded size = 3 and
      # silently ignored its own point_size argument. All call sites pass the
      # default (3), so this is a behavioral no-op in practice.
      geom_point(shape = 21, size = point_size, fill = "black", alpha = 0.7) +
      geom_abline(slope = 1, color = "red")
  } else {
    cor_plot <- cor_plot +
      geom_hex(aes(fill = after_stat(log10(count)))) +
      scale_fill_viridis_c() +
      geom_abline(slope = 1, color = "red")
  }

  detection_plot <- TPM_df %>%
    filter((mean_TPM_AoU > detection_threshold & mean_TPM_GTEx < detection_threshold) |
           (mean_TPM_AoU < detection_threshold & mean_TPM_GTEx > detection_threshold)) %>%
    mutate(type = case_when(mean_TPM_AoU > mean_TPM_GTEx ~ "AoU", TRUE ~ "GTEx")) %>%
    dplyr::count(type) %>%
    ggplot(aes(x = type, y = n, fill = type)) +
    geom_col(show.legend = FALSE) +
    ylab(axis_label) + xlab("Cohort") +
    scale_fill_carto_d(palette = "Safe")

  cor_plot + detection_plot
}

# =============================================================================
# Builders
# =============================================================================
gene_annotations <- function(refresh = FALSE)
  build_or_cache("gene_annotations", refresh = refresh, {
    rtracklayer::import(input_path("gencode_gtf")) %>%
      data.frame() %>% select(gene_id, gene_type) %>% distinct()
  })

# Pure join of two compute_mean_TPM() outputs: version-strip the gene id,
# inner-join AoU to GTEx, set missing means to 0. Factored out so it is testable
# without the GTF import / GCT reads that merged_tpms() performs.
join_mean_tpms <- function(aou_means, gtex_means) {
  strip <- function(d) d %>% distinct() %>% mutate(gene_id = str_remove(Name, "\\..*"))
  strip(aou_means) %>%
    select(-Name) %>%
    inner_join(strip(gtex_means) %>% select(-Description, -gene_type, -Name),
               by = "gene_id", suffix = c("_AoU", "_GTEx")) %>%
    mutate(mean_TPM_GTEx = coalesce(mean_TPM_GTEx, 0),
           mean_TPM_AoU  = coalesce(mean_TPM_AoU, 0))
}

# AoU + GTEx per-gene mean TPM, joined and gene-type-attached. The GCT files have
# a 2-line header, hence skip = 2.
merged_tpms <- function(refresh = FALSE)
  build_or_cache("merged_tpms", refresh = refresh, {
    ann  <- gene_annotations()
    aou  <- as_tibble(data.table::fread(input_path("aou_tpms"),  skip = 2)) %>% compute_mean_TPM(ann)
    gtex <- as_tibble(data.table::fread(input_path("gtex_tpms"), skip = 2)) %>% compute_mean_TPM(ann)
    join_mean_tpms(aou, gtex)
  })
