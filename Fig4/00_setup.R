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
