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
