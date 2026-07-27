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

# =============================================================================
# Panels A-C (population / constraint view). Each build_panel*() returns a single
# UNTAGGED object; tagging is applied in fig5/assemble.R.
#   A  fine-mapped rare ClinGen variants by ancestry (AD/AR) + ClinVar proportion
#   B  |average QTL effect| by variant-effect category (Control/AD/AR)
#   C  log2 aFC of the strongest decreasing variant per gene (AD/AR)
# =============================================================================

# Rare, confidently fine-mapped variants at ClinGen vs low-shet control genes. The
# gene_assignments left_join fans a gene curated under BOTH AD and AR into two rows
# (counted under both classes) — preserved from the original; see NOTES.md.
variant_effect_by_inheritance <- function(eqtl, clingen) {
  gene_assignments <- clingen %>%
    distinct(targetFromSourceId, allelicRequirements) %>%
    filter(allelicRequirements != "SD")
  eqtl %>%
    mutate(ControlGene = post_mean < 0.01 & ClinGen == FALSE) %>%
    filter((ClinGen == TRUE | ControlGene == TRUE), pip > 0.9, rare == TRUE) %>%
    select(variant, log2_aFC, gene_name, posterior_mean) %>%
    left_join(gene_assignments, by = c("gene_name" = "targetFromSourceId")) %>%
    distinct() %>%
    mutate(allelicRequirements = factor(replace_na(allelicRequirements, "Control"),
                                        levels = c("Control", "AD", "AR")))
}

build_panelB <- function(eqtl, clingen) {
  vebi <- variant_effect_by_inheritance(eqtl, clingen)
  sign_summary <- vebi %>%
    mutate(dir = if_else(log2_aFC < 0, "Decreasing", "Increasing")) %>%
    group_by(allelicRequirements, dir) %>%
    summarise(MeanEffect = mean(abs(posterior_mean), na.rm = TRUE),
              se = sd(abs(posterior_mean), na.rm = TRUE) / sqrt(n()), .groups = "drop")
  aggregate_summary <- vebi %>%
    group_by(allelicRequirements) %>%
    summarise(MeanEffect = mean(abs(posterior_mean), na.rm = TRUE),
              se = sd(abs(posterior_mean), na.rm = TRUE) / sqrt(n()), .groups = "drop") %>%
    mutate(dir = "Aggregate")
  average_effect_summary <- bind_rows(sign_summary, aggregate_summary) %>%
    mutate(dir = factor(dir, levels = c("Aggregate", "Decreasing", "Increasing")),
           conf.low = MeanEffect - 1.96 * se, conf.high = MeanEffect + 1.96 * se)
  category_lines <- tibble(x = c(0.60, 1.60, 2.60), xend = c(1.40, 2.40, 3.40))

  ggplot(average_effect_summary, aes(dir, MeanEffect, color = allelicRequirements)) +
    geom_segment(data = category_lines, aes(x = x, xend = xend, y = -Inf, yend = -Inf),
                 inherit.aes = FALSE, color = "grey72", linewidth = 0.55, lineend = "round") +
    geom_linerange(aes(ymin = conf.low, ymax = conf.high),
                   position = position_dodge(width = 0.5), linewidth = 0.6) +
    geom_point(position = position_dodge(width = 0.5), size = 2.3) +
    scale_x_discrete(drop = FALSE) +
    scale_color_manual(values = inheritance_colors, drop = FALSE) +
    labs(x = "Variant Effect Category", y = "|Average QTL effect|", color = NULL) +
    coord_cartesian(clip = "off") +
    theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.x = element_text(margin = margin(t = 1)),
          legend.position = c(0.5, 0.96), legend.direction = "horizontal")
}

build_panelC <- function(eqtl, clingen, nohom = NULL, n_genes = 7,
                         required_genes = c("HEXA", "MLC1")) {
  vebi <- variant_effect_by_inheritance(eqtl, clingen)
  candidate <- vebi %>%
    {if (!is.null(nohom)) filter(., variant %in% nohom$variant) else .} %>%
    filter(allelicRequirements %in% c("AD", "AR"), !is.na(log2_aFC), log2_aFC < 0) %>%
    group_by(gene_name, allelicRequirements) %>%
    slice_min(log2_aFC, n = 1, with_ties = FALSE) %>% ungroup() %>%
    mutate(required_gene = gene_name %in% required_genes)
  high_effect <- candidate %>%
    group_by(allelicRequirements) %>%
    arrange(desc(required_gene), log2_aFC, .by_group = TRUE) %>%
    slice_head(n = n_genes) %>% ungroup() %>%
    mutate(allelicRequirements = factor(allelicRequirements, levels = c("AD", "AR")),
           display_gene = fct_reorder(gene_name, log2_aFC, .desc = TRUE))

  ggplot(high_effect, aes(log2_aFC, display_gene, color = allelicRequirements)) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.5) +
    geom_segment(aes(x = 0, xend = log2_aFC, y = display_gene, yend = display_gene),
                 color = "grey78", linewidth = 0.65) +
    geom_point(size = 3.2) +
    facet_grid(allelicRequirements ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_y_discrete(position = "right") +
    scale_color_manual(values = inheritance_colors, breaks = c("AD", "AR"), guide = "none") +
    scale_x_continuous(limits = c(-7, 0), breaks = c(-6, -4, -2, 0),
                       expand = expansion(mult = c(0.03, 0.02))) +
    labs(x = expression("log"[2]*" aFC (decreasing variants only)"), y = NULL) +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(face = "bold", angle = 0),
          axis.text.y.right = element_text(hjust = 0),
          panel.spacing.y = unit(1, "lines"),
          axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          plot.margin = margin(5.5, 20, 5.5, 5.5))
}

build_panelA <- function(eqtl, clingen) {
  clingen_only <- eqtl %>%
    filter(ClinGen == TRUE) %>%
    left_join(clingen, by = c("gene_id_nov" = "targetId")) %>%
    mutate(causal = !is.na(log2_aFC))
  rare_clingen <- clingen_only %>%
    filter(rare == TRUE) %>% mutate(InClinvar = !is.na(clinvar_classification))

  ancestry_count <- rare_clingen %>%
    count(allelicRequirements, gvs_max_subpop, InClinvar) %>%
    mutate(gvs_max_subpop = str_remove_all(gvs_max_subpop, "gvs_|_af"),
           gvs_max_subpop = if_else(gvs_max_subpop == "oth", "admix", gvs_max_subpop)) %>%
    filter(allelicRequirements != "SD", !is.na(gvs_max_subpop))

  ancestry_plot <- ancestry_count %>%
    ggplot(aes(gvs_max_subpop, n, fill = interaction(gvs_max_subpop, InClinvar))) +
    geom_col(color = "transparent") +
    geom_point(aes(color = InClinvar), x = Inf, y = Inf, inherit.aes = FALSE) +
    facet_wrap(~allelicRequirements) +
    scale_fill_manual(values = ancestry_fill, guide = "none") +
    scale_x_discrete(labels = toupper) +
    scale_color_manual(values = c(`FALSE` = "grey70", `TRUE` = "grey30"),
                       labels = c(`FALSE` = "Novel Variant", `TRUE` = "ClinVar Annotated"),
                       name = NULL) +
    guides(color = guide_legend(override.aes = list(size = 5, shape = 15))) +
    geom_segment(aes(x = 0.6, xend = 7.4, y = Inf, yend = Inf), inherit.aes = FALSE,
                 color = "grey70", linewidth = 0.6) +
    labs(y = "Fine-mapped Rare Variants") +
    theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title.x = element_blank(), legend.position = c(0.9, 0.85))

  clinvar_prop <- rare_clingen %>%
    mutate(clinvar_class = case_when(
      str_detect(clinvar_classification, "pathogenic") ~ "P/LP",
      str_detect(clinvar_classification, "uncertain significance") ~ "VUS",
      str_detect(clinvar_classification, "benign") ~ "B/LB",
      TRUE ~ "Unannotated")) %>%
    count(allelicRequirements, clinvar_class) %>%
    filter(allelicRequirements != "SD") %>%
    group_by(allelicRequirements) %>% mutate(prop = n / sum(n)) %>% ungroup() %>%
    ggplot(aes(prop, allelicRequirements, fill = clinvar_class)) +
    geom_col(position = "stack", width = 0.85) +
    facet_wrap(~allelicRequirements, nrow = 1, scales = "free_y") +
    paletteer::scale_fill_paletteer_d("futurevisions::mars") +
    scale_x_continuous(breaks = c(0, .25, .5, .75, 1),
                       labels = \(x) paste0(x * 100, "%"), expand = expansion(mult = 0)) +
    scale_y_discrete(expand = expansion(add = 0.06)) +
    guides(fill = guide_legend(nrow = 1)) +
    labs(fill = NULL, x = NULL) +
    theme(axis.text.y = element_blank(), axis.line.y = element_blank(),
          axis.title.y = element_blank(), axis.ticks.y = element_blank(),
          axis.line.x = element_blank(),
          axis.text.x = element_text(size = 9.5, margin = margin(t = 0)),
          axis.ticks.x = element_blank(),
          strip.text = element_blank(), panel.spacing.x = unit(1.6, "lines"),
          plot.margin = margin(0, 8, 0, 8),
          legend.position = "bottom", legend.direction = "horizontal",
          legend.box.spacing = unit(2, "pt"), legend.margin = margin(t = 0, b = 0),
          legend.text = element_text(size = 9), legend.title = element_text(size = 9),
          legend.key.width = unit(0.45, "cm"), legend.key.height = unit(0.35, "cm"))

  ancestry_plot / clinvar_prop + plot_layout(heights = c(7, 1))
}
