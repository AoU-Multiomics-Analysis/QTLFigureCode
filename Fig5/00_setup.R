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

# =============================================================================
# Data loaders (cached). GeneBayes shet + ClinGen Definitive/Strong genes, then
# the eQTL SuSiE table annotated with shet + ClinGen flag + rarity.
#
# Each heavy loader wraps a PURE seam (frames -> frames) so the join / filter
# logic is unit-testable without the real Excel / gzip reads.
# =============================================================================

# ClinGen Definitive/Strong AD/AR/SD genes mapped to Ensembl ids via the SuSiE
# gene map. NB: a gene curated under BOTH AD and AR yields two rows here; that
# fan-out is deliberately preserved downstream (Panels B/C) — see NOTES.md.
derive_clingen <- function(clingen_raw, susie_gene_map) {
  clingen_raw %>%
    filter(confidence %in% c("Definitive", "Strong"), moi %in% c("AD", "AR", "SD")) %>%
    distinct(gene, moi, confidence) %>%
    transmute(targetFromSourceId = gene, allelicRequirements = moi, confidence) %>%
    left_join(susie_gene_map, by = c("targetFromSourceId" = "gene_name")) %>%
    dplyr::rename(targetId = gene_id)
}

load_reference <- function() {
  GeneBayes <- read_excel(input_path("genebayes"), sheet = "Supplementary Table 1") %>%
    transmute(ensg, post_mean)
  clingen_raw <- fread(input_path("clingen"), skip = 6, header = FALSE, select = 1:7,
    col.names = c("gene", "hgnc", "disease", "mondo", "moi", "sop", "confidence"))
  susie_gene_map <- fread(input_path("susie"), select = c("gene_id", "gene_name")) %>%
    distinct() %>% mutate(gene_id = sub("\\..*", "", gene_id))
  list(GeneBayes = GeneBayes, clingen = derive_clingen(clingen_raw, susie_gene_map))
}

# Annotate the SuSiE eQTL table with GeneBayes shet, a ClinGen flag, and rarity.
# `clingen_genes` is passed as a plain character vector (not reference$...) so it
# cannot collide with the SuSiE table's own `ref` column under data masking.
annotate_eqtl <- function(susie, genebayes, clingen_genes) {
  susie %>%
    mutate(gene_id_nov = sub("\\..*", "", gene_id)) %>%
    left_join(genebayes, by = c("gene_id_nov" = "ensg")) %>%
    mutate(ClinGen = gene_name %in% clingen_genes, rare = gvs_max_af < 0.01)
}

load_eqtl_susie <- function(reference)
  annotate_eqtl(fread(input_path("susie")),
                reference$GeneBayes, reference$clingen$targetFromSourceId)

build_variant_freq <- function(eqtl) {
  eqtl %>%
    transmute(variant_key = normalize_variant_id(variant),
              max_af = gvs_max_af, subpop = gvs_max_subpop, rare, log2_aFC) %>%
    distinct(variant_key, .keep_all = TRUE) %>%
    mutate(freq_label = make_freq_label(max_af, subpop, rare),
           afc_label  = sprintf("log₂ aFC %+.1f", log2_aFC))
}

# optional private list refining Panel C; NULL (with a note) if absent
load_nohomozygotes <- function() {
  f <- tryCatch(input_path("nohom"), error = function(e) NA_character_)
  if (!is.na(f) && file.exists(f))
    return(fread(f) %>% mutate(variant = paste(V1, V2, V3, V4, sep = "_")))
  message("NOTE: nohom input not found (private) — Panel C ranks the strongest ",
          "decreasing variant per gene without it.")
  NULL
}

# cached getters used by the panel scripts
get_reference    <- function(refresh = FALSE) build_or_cache("reference", load_reference(), refresh)
get_eqtl         <- function(refresh = FALSE) build_or_cache("eqtl_susie", load_eqtl_susie(get_reference()), refresh)
get_variant_freq <- function(refresh = FALSE) build_or_cache("variant_freq", build_variant_freq(get_eqtl()), refresh)
