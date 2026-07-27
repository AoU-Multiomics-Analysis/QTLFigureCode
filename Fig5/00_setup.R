# =============================================================================
# Figure 5 — shared setup (sourced by every panel / assembly script).
# Run all scripts from the PROJECT ROOT. Each script starts with:
#   source("fig5/00_setup.R")
#
# Figure 5 covers rare-variant eQTL constraint and sequence-model variant-effect
# dissection: fine-mapped rare ClinGen variants (A-C), and two worked examples
# (ACVRL1 promoter, PromoterAI: D/E; SPR element, Borzoi browser: F/G).
# =============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(readxl)
  library(patchwork)
  library(paletteer)
  library(rcartocolor)
  library(colorspace)
  library(scales)
  library(ggrepel)
  library(cowplot)
  library(ggseqlogo)
  library(showtext)
  library(sysfonts)
})
# rtracklayer is used namespaced (rtracklayer::import) so it is NOT attached —
# attaching it lets S4Vectors/BiocGenerics mask dplyr verbs (rename/slice/desc).

# --- Paths (override via env vars if needed) --------------------------------
INPUT    <- Sys.getenv("FIG5_INPUT",    ".")             # staging root (repo root)
DATA_DIR <- Sys.getenv("FIG5_DATA",     "fig5/data")     # memoized derived tables
OBJ_DIR  <- Sys.getenv("FIG5_OBJ",      "fig5/objects")  # per-panel ggplot objects
OUTDIR   <- Sys.getenv("FIG5_OUTDIR",   ".")             # final Figure5.pdf/png
MANIFEST <- Sys.getenv("FIG5_MANIFEST", "fig5/inputs.tsv")
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
SZ_AXIS_TITLE <- 12; SZ_AXIS_TEXT <- 10; SZ_STRIP <- 11; SZ_TAG <- 20; SZ_TITLE <- 12

theme_set(theme_classic())
theme_update(
  text             = element_text(family = FONT),
  axis.title       = element_text(size = SZ_AXIS_TITLE),
  axis.text        = element_text(size = SZ_AXIS_TEXT),
  strip.text       = element_text(size = SZ_STRIP),
  strip.background = element_blank(),
  plot.tag         = element_text(size = SZ_TAG, face = "bold"),
  plot.margin      = margin(9, 12, 9, 12))

# --- Palettes ---------------------------------------------------------------
# inheritance / gene class (Panels A-C)
inheritance_colors <- c(Control = "#D43F3A", AD = "#EEA236", AR = "#5CB85C")

# ancestry base colours (Panel A), named by sub-population
ancestry_base <- setNames(carto_pal(7, "Temps"),
                          c("afr", "amr", "eas", "eur", "mid", "admix", "sas"))
# Panel A stacked-bar fill: lighter = novel, darker = ClinVar-annotated
ancestry_fill <- unlist(lapply(names(ancestry_base), function(pop)
  c(setNames(lighten(unname(ancestry_base[pop]), 0.15), paste0(pop, ".FALSE")),
    setNames(darken(unname(ancestry_base[pop]), 0.10),  paste0(pop, ".TRUE")))))

# ENCODE cCRE shading (Panels D/F): PLS = promoter-like, ELS = enhancer-like
CCRE_FILL <- c(PLS = "#E15759", ELS = "#F4A261")
# alternate-allele nucleotide colours (Panel E, PromoterAI)
nucleotide_colors <- c(A = "#E64B35", C = "#00A087", G = "#F39B7F", T = "#4DBBD5")

# credible-set palette: a named Prism vector for a set of cs_id levels
cs_palette <- function(levels) {
  levels <- unique(levels[!is.na(levels)])
  setNames(carto_pal(max(3, length(levels)), "Prism")[seq_along(levels)], levels)
}

# --- Small utilities --------------------------------------------------------
# signed-square y-axis transform (Panel E): compresses tails so small
# base-resolution effects stay visible next to the causal spike
signed_square_trans <- function()
  scales::trans_new("signed_square",
    transform = function(x) sign(x) * x^2,
    inverse   = function(x) sign(x) * sqrt(abs(x)))

fmt_af <- function(af) dplyr::case_when(
  is.na(af)   ~ "NA",
  af >= 0.01  ~ sprintf("%.1f%%", af * 100),
  af >= 0.001 ~ sprintf("%.2f%%", af * 100),
  TRUE        ~ sprintf("%.3f%%", af * 100))

make_freq_label <- function(max_af, subpop, rare) {
  pop  <- ifelse(is.na(subpop), "", sprintf(" (%s)", toupper(subpop)))
  base <- sprintf("MAF %s%s", fmt_af(max_af), pop)
  ifelse(rare, paste0("rare · ", base), base)
}

normalize_variant_id <- function(x)
  str_replace_all(str_replace_all(as.character(x), "[:*]", "_"), "_+", "_")
make_variant_id <- function(chrom, pos, ref, alt)
  normalize_variant_id(paste(chrom, pos, toupper(ref), toupper(alt), sep = "_"))

# --- Derived-data cache (identical contract to Fig4) ------------------------
# Lazy `expr`, cache to DATA_DIR/<name>.rds, keyed by NAME only. Clear DATA_DIR
# after changing an input.
build_or_cache <- function(name, expr, refresh = FALSE) {
  f <- file.path(DATA_DIR, paste0(name, ".rds"))
  if (!refresh && file.exists(f)) return(readRDS(f))
  message("building ", name, " ...")
  val <- force(expr)
  saveRDS(val, f)
  val
}
