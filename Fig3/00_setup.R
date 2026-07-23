# =============================================================================
# Figure 3 — shared setup (sourced by every panel / supplemental / assembly).
# Run all scripts from the PROJECT ROOT. Each script starts with:
#   source("fig3/00_setup.R")
# Re-sourcing is cheap and idempotent.
#
# Figure 3 covers cross-QTL colocalization (FastENLOC) between each pair of
# molecular assays — eQTL-pQTL, sQTL-pQTL, eQTL-sQTL — at gene, credible-set
# (signal) and variant level.
#
# See docs/superpowers/specs/2026-07-23-fig3-cleanup-design.md
# =============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(patchwork)
  library(igraph)
  library(broom)
  library(scales)
  library(showtext)
  library(sysfonts)
  library(paletteer)
})

# --- Paths (override via env vars if needed) --------------------------------
# INPUT is the staging root; it defaults to the current working directory, so
# gsutil-staged files land under the cwd in their manifest subpaths (enloc/,
# finemapping/, ...). Set FIG3_INPUT to keep inputs in a separate folder.
INPUT    <- Sys.getenv("FIG3_INPUT",    ".")            # staging root (cwd)
DATA_DIR <- Sys.getenv("FIG3_DATA",     "fig3/data")    # memoized derived tables
OBJ_DIR  <- Sys.getenv("FIG3_OBJ",      "fig3/objects") # per-panel ggplot objects
OUTDIR   <- Sys.getenv("FIG3_OUTDIR",   ".")            # final Figure3.pdf/png
MANIFEST <- Sys.getenv("FIG3_MANIFEST", "fig3/inputs.tsv")
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
# Figure 3 is an 18x12" canvas, so these are much smaller than Figure 1's.
SZ_AXIS_TITLE <- 11; SZ_AXIS_TEXT <- 9;  SZ_STRIP  <- 12
SZ_LEG_TEXT   <- 9;  SZ_LEG_TITLE <- 10; SZ_TAG    <- 14
SZ_TITLE      <- 11
SZ_ANNOT      <- 2.7           # mm — geom_label()/annotate() only

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

# The three assay pairs, in the order used throughout the figure.
COLOC_PAIRS <- c("eQTL-pQTL", "sQTL-pQTL", "eQTL-sQTL")

# Panel A fill — number of colocalized molecular layers.
LAYER_FILL <- c("0" = "#D9D9D9", "1" = "#A9C1D9", "2" = "#7B9AB8", "3" = "#577590")

# Panel B fill — number of colocalized signals.
SIGNAL_FILL <- c("0" = "#E9DCDC", "1" = "#D8BBBB", "2" = "#C89A9A", "3" = "#B37979",
                 "4" = "#9E5E5E", "5" = "#8A4D4D", "6" = "#763F3F", "7" = "#643232")

# Per-ancestry allele-frequency columns of the AoU annotation table, used to
# derive a per-variant maximum MAF. Named explicitly so a schema change fails
# loudly rather than silently changing which columns are folded.
AF_POP_COLS <- c("gvs_afr_af", "gvs_amr_af", "gvs_eas_af", "gvs_eur_af",
                 "gvs_mid_af", "gvs_sas_af", "gvs_oth_af")

# =============================================================================
# FastENLOC output readers
# -----------------------------------------------------------------------------
# The .sig.out and .snp.out files are whitespace-aligned and their header line
# does not line up with the data columns, so they are parsed positionally: drop
# the header, split each row on runs of whitespace, and apply known names.
# (.gene.out is a normal TSV and is read with fread.)
# =============================================================================
ENLOC_CS_COLS      <- c("Signal", "Num_SNP", "CPIP_qtl", "CPIP_gwas_marginal",
                        "CPIP_gwas_qtl_prior", "RCP", "LCP")
ENLOC_VARIANT_COLS <- c("Signal", "SNP", "PIP_qtl", "PIP_gwas_marginal",
                        "PIP_gwas_qtl_prior", "SCP")

LoadEnlocTable <- function(path, cols) {
  x <- readLines(path)
  x <- x[-1]                                   # drop the header line
  x <- x[nzchar(trimws(x))]                    # drop blank lines
  mat <- do.call(rbind, strsplit(trimws(x), "\\s+"))
  if (ncol(mat) != length(cols))
    stop("Expected ", length(cols), " columns in ", path, ", found ", ncol(mat))
  out <- as.data.frame(mat, stringsAsFactors = FALSE)
  names(out) <- cols
  # strsplit() yields character; convert here so downstream code needs no
  # scattered as.numeric() / across(contains('PIP'), as.numeric) calls.
  tibble::as_tibble(type.convert(out, as.is = TRUE))
}

# Label rows as colocalized by controlling the Bayesian FDR at `control_level`:
# sort local FDRs (1 - posterior) ascending, take the largest prefix whose
# running mean stays at or below the target, and threshold at its posterior.
#
# `prob_col` and `label_col` are plain COLUMN-NAME STRINGS. The original used
# rlang::ensym() here, but every call is internal and passes a name held in a
# variable, which is exactly the case non-standard evaluation cannot handle.
#
# The original stashed the threshold in attr(out, "threshold"), which the
# following dplyr verbs silently dropped; it is reported via message() instead.
add_fdr_threshold_label <- function(df, prob_col, control_level = 0.05,
                                    label_col = "Coloc", quiet = FALSE) {
  stopifnot(is.character(prob_col), length(prob_col) == 1L)
  if (!prob_col %in% names(df))
    stop("Column '", prob_col, "' not found; have: ", paste(names(df), collapse = ", "))

  probs <- as.numeric(df[[prob_col]])

  lfdr  <- sort(1 - probs, na.last = NA)
  FDR   <- cumsum(lfdr) / seq_along(lfdr)
  # Tolerance guards an exact-boundary case (e.g. probs = 0.95, control_level =
  # 0.05) where `1 - 0.95` rounds to a double a hair above 0.05, which would
  # otherwise make a should-pass row fail `FDR <= control_level`.
  valid <- which(FDR <= control_level + 1e-8)

  thresh <- if (length(valid) == 0) Inf else 1 - lfdr[max(valid)]

  if (!quiet)
    message(sprintf("    %s >= %.4g at FDR %.0f%%  (%d / %d rows labelled)",
                    prob_col, thresh, 100 * control_level,
                    sum(probs >= thresh, na.rm = TRUE), length(probs)))

  df[[label_col]] <- probs >= thresh
  df
}
