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
OUTDIR   <- Sys.getenv("FIG5_OUTDIR", "fig5")         # final Figure5.pdf/png
MANIFEST <- Sys.getenv("FIG5_MANIFEST", "fig5/inputs.tsv")
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OBJ_DIR,  showWarnings = FALSE, recursive = TRUE)

dir.create(OUTDIR,   showWarnings = FALSE, recursive = TRUE)

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

# =============================================================================
# Shared fine-mapping plotting (Panels D & F): cCRE bands, gene-model track, and
# the regional-PIP + zoom-inset plotter. Returned objects are UNTAGGED; tags are
# applied in fig5/assemble.R.
# =============================================================================
.gene_exons_all <- fread(input_path("gene_exons"), col.names = c("gene", "chrom", "start", "end", "strand"))
.gene_spans_all <- fread(input_path("gene_spans"), col.names = c("gene", "chrom", "start", "end", "strand"))

# ENCODE cCRE shading: keep only cCREs overlapping a fine-mapped variant. Light
# fill + crisp vertical edges only, so adjacent elements stay distinct.
ccre_bands <- function(ccre, scale = 1e6, alpha = 0.16, vars = NULL) {
  if (is.null(ccre) || nrow(ccre) == 0) return(NULL)
  if (!is.null(vars))
    ccre <- ccre |> filter(map2_lgl(start, end, \(s, e) any(vars >= s & vars <= e)))
  if (nrow(ccre) == 0) return(NULL)
  ccre  <- ccre |>
    mutate(fillc = if_else(grepl("PLS", type), CCRE_FILL["PLS"], CCRE_FILL["ELS"])) |>
    select(start, end, fillc)                 # drop `type` (collides with facet_wrap(~type))
  edges <- ccre |> pivot_longer(c(start, end), values_to = "x")
  list(
    geom_rect(data = ccre, aes(xmin = start / scale, xmax = end / scale),
              ymin = -Inf, ymax = Inf, fill = ccre$fillc, alpha = alpha, inherit.aes = FALSE),
    geom_segment(data = edges, aes(x = x / scale, xend = x / scale), y = -Inf, yend = Inf,
                 colour = edges$fillc, linewidth = 0.3, inherit.aes = FALSE))
}

# gene-model track sharing x limits with the PIP panel above it
gene_model_track <- function(xlim_mb, scale = 1e6, exon_h = 0.34, label_size = 3.4,
                             exons = .gene_exons_all, spans = .gene_spans_all,
                             expand = expansion(mult = 0), x_title = "Genomic position (Mb)",
                             x_labeller = waiver(), right_strip = TRUE) {
  win <- xlim_mb * scale
  sp  <- spans |> filter(start < win[2], end > win[1])
  ex  <- exons |> filter(gene %in% sp$gene, start < win[2], end > win[1])
  lab <- sp |> mutate(tss = if_else(strand == "+", start, end),
                      arr = if_else(strand == "+", 1, -1) * diff(win) * 0.045)
  p <- ggplot() +
    geom_segment(data = sp, aes(x = pmax(start / scale, xlim_mb[1]),
                                xend = pmin(end / scale, xlim_mb[2]), y = 0, yend = 0),
                 colour = "grey45", linewidth = 0.5) +
    geom_rect(data = ex, aes(xmin = start / scale, xmax = end / scale,
                             ymin = -exon_h, ymax = exon_h), fill = "grey35") +
    geom_segment(data = lab, aes(x = tss / scale, xend = tss / scale,
                                 y = 0, yend = exon_h * 2.2), colour = "grey45", linewidth = 0.5) +
    geom_segment(data = lab, aes(x = tss / scale, xend = (tss + arr) / scale,
                                 y = exon_h * 2.2, yend = exon_h * 2.2), colour = "grey45",
                 linewidth = 0.5, arrow = arrow(length = unit(3.5, "pt"), type = "closed")) +
    geom_text(data = lab, aes(x = tss / scale, y = 0, label = gene),
              size = label_size, fontface = "italic", hjust = 1,
              nudge_x = -diff(xlim_mb) * 0.012) +
    scale_x_continuous(limits = xlim_mb, expand = expand, labels = x_labeller) +
    scale_y_continuous(limits = c(-exon_h * 1.6, exon_h * 2.9)) +
    coord_cartesian(clip = "off") +
    labs(x = x_title, y = NULL) +
    theme(panel.background = element_blank(), panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.x = element_text(size = 11), axis.title.x = element_text(size = 12),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.line.y = element_blank(), axis.title.y = element_blank(),
          legend.position = "none", plot.margin = margin(1, 8, 4, 8))
  if (right_strip)
    p <- p + facet_wrap(~"x", strip.position = "right") +
      theme(strip.text = element_text(size = 12, colour = NA),
            strip.background = element_blank(), panel.spacing = unit(0, "pt"))
  p
}

# regional PIP plot with a zoom inset and an aligned gene-model track (untagged)
plot_finemapping_with_inset <- function(finemapping_df, pip_thresh = 0.9,
    zoom_center, zoom_window_bp = 500,
    inset_left = 0.24, inset_bottom = 0.05, inset_right = 0.96, inset_top = 0.76,
    box_padding_bp = 1000, box_padding_pip = 0.05, cs_colors,
    panel_title = NULL, ccre = NULL, x_limits = NULL) {

  plot_df <- finemapping_df %>%
    filter(!is.na(position), !is.na(pip)) %>%
    mutate(causal = pip >= pip_thresh, cs_id = factor(cs_id),
           in_zoom = abs(position - zoom_center) <= zoom_window_bp)
  zoom_df   <- plot_df %>% filter(in_zoom)
  iz_causal <- plot_df$position[plot_df$causal &  plot_df$in_zoom]
  oz_causal <- plot_df$position[plot_df$causal & !plot_df$in_zoom]

  regional_x_range  <- range(plot_df$position / 1e6, na.rm = TRUE)
  regional_x_pad    <- max(diff(regional_x_range) * 0.02, 0.001)
  regional_x_limits <- if (is.null(x_limits)) regional_x_range + c(-1, 1) * regional_x_pad else x_limits
  regional_y_limits <- c(0, 1.08)

  boxed <- zoom_df %>% filter(causal); if (nrow(boxed) == 0) boxed <- zoom_df
  box_start_mb  <- max(min(boxed$position) - box_padding_bp, min(plot_df$position)) / 1e6
  box_end_mb    <- min(max(boxed$position) + box_padding_bp, max(plot_df$position)) / 1e6
  box_center_mb <- mean(c(box_start_mb, box_end_mb))
  box_y_min <- max(regional_y_limits[1], min(boxed$pip) - box_padding_pip)
  box_y_max <- min(regional_y_limits[2], max(boxed$pip) + box_padding_pip)
  inset_top_data <- regional_y_limits[1] + inset_top * diff(regional_y_limits)

  regional_plot <- ggplot(plot_df, aes(position / 1e6, pip)) +
    ccre_bands(ccre, scale = 1e6, vars = oz_causal) +
    geom_hline(yintercept = pip_thresh, linetype = "dashed", linewidth = 0.5) +
    annotate("segment", x = box_center_mb, xend = box_center_mb,
             y = box_y_min, yend = inset_top_data, color = "grey55", linewidth = 0.4) +
    annotate("rect", xmin = box_start_mb, xmax = box_end_mb, ymin = box_y_min, ymax = box_y_max,
             fill = "grey90", alpha = 0.18, color = "grey35", linewidth = 0.5) +
    geom_point(data = filter(plot_df, !causal), aes(fill = cs_id), shape = 21,
               size = 1.8, stroke = 0.35, color = "grey55", alpha = 0.5, show.legend = FALSE) +
    geom_point(data = filter(plot_df, causal), aes(fill = cs_id), shape = 23,
               size = 3.2, stroke = 0.65, color = "black", show.legend = FALSE) +
    geom_point(data = filter(plot_df, causal & rare & !in_zoom), shape = 5, size = 4.8,
               colour = "#8C6D1F", stroke = 0.9, show.legend = FALSE) +
    geom_point(data = tibble(type = plot_df$type[1]),
               aes(x = regional_x_limits[1], y = -0.5, shape = "rare (MAF < 1%)"),
               colour = "#8C6D1F", size = 3.2, stroke = 0.9, inherit.aes = FALSE) +
    scale_shape_manual(name = NULL, values = c("rare (MAF < 1%)" = 5)) +
    facet_wrap(~type, ncol = 1, strip.position = "right") +
    scale_fill_manual(values = cs_colors, na.value = "grey80") +
    scale_x_continuous(limits = regional_x_limits, expand = expansion(mult = 0)) +
    scale_y_continuous(breaks = seq(0, 1, 0.25), limits = regional_y_limits,
                       expand = expansion(mult = 0)) +
    labs(x = NULL, y = "PIP", title = panel_title) +
    theme(strip.text = element_text(size = 12),
          panel.border = element_rect(fill = NA, color = "black", linewidth = 0.6),
          axis.line = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "top", legend.justification = "left",
          legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(1, "pt"),
          legend.text = element_text(size = 9), legend.key.size = unit(4, "mm"),
          legend.key = element_blank(),
          plot.title = element_text(size = 12.5, hjust = 0.5, face = "bold", margin = margin(b = 4)),
          plot.margin = margin(6, 8, 2, 8))

  zoom_plot <- ggplot(zoom_df, aes(position, pip)) +
    ccre_bands(ccre, scale = 1, vars = iz_causal) +
    geom_hline(yintercept = pip_thresh, linetype = "dashed", linewidth = 0.4) +
    geom_point(data = filter(zoom_df, !causal), aes(fill = cs_id), shape = 21,
               size = 1.7, stroke = 0.3, color = "grey55", alpha = 0.5, show.legend = FALSE) +
    geom_point(data = filter(zoom_df, causal), aes(fill = cs_id), shape = 23,
               size = 3, stroke = 0.6, color = "black", show.legend = FALSE) +
    geom_point(data = filter(zoom_df, causal & rare), shape = 5, size = 4.6,
               colour = "#8C6D1F", stroke = 0.85, show.legend = FALSE) +
    scale_fill_manual(values = cs_colors, na.value = "grey80") +
    scale_x_continuous(
      breaks = (\(b) b[seq(1, length(b), 2)])(pretty(zoom_center + c(-1, 1) * zoom_window_bp, 4)),
      limits = c(zoom_center - zoom_window_bp, zoom_center + zoom_window_bp),
      labels = \(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE),
      expand = expansion(mult = c(0.04, 0.04))) +
    scale_y_continuous(breaks = c(0, 0.5, pip_thresh), limits = c(0, 1.12),
                       expand = expansion(mult = c(0.02, 0.08))) +
    labs(x = NULL, y = NULL) +
    theme(panel.border = element_rect(fill = NA, color = "grey35", linewidth = 0.5),
          axis.line = element_blank(), axis.text = element_text(size = 8.5),
          legend.position = "none",
          plot.background = element_rect(fill = "white", color = NA))

  main   <- regional_plot +
    inset_element(zoom_plot, left = inset_left, bottom = inset_bottom,
                  right = inset_right, top = inset_top, align_to = "panel",
                  on_top = TRUE, ignore_tag = TRUE)
  gtrack <- gene_model_track(regional_x_limits, expand = expansion(mult = 0))
  main / gtrack + plot_layout(heights = c(1, 0.17))
}

# =============================================================================
# Panels D/E (ACVRL1 worked example). Untagged; tags applied in assemble.R.
# =============================================================================
ACVRL1_FOCAL <- 51912316   # ACVRL1 promoter variant (chr12, hg38)

# PromoterAI for one transcript around a focal position, joined to the QTL
prepare_promoterAI_locus <- function(promoter_data, qtl_data, gene_symbol, transcript,
                                     focal_position, window_bp = 100) {
  qtl_lookup <- qtl_data %>%
    transmute(variant_key = normalize_variant_id(variant), pip, causal, cs_id) %>%
    distinct(variant_key, .keep_all = TRUE)
  promoter_data %>%
    filter(gene == gene_symbol, transcript_id == transcript) %>%
    mutate(relative_position = pos - focal_position,
           variant     = make_variant_id(chrom, pos, ref, alt),
           variant_key = normalize_variant_id(variant)) %>%
    filter(abs(relative_position) < window_bp) %>%
    left_join(qtl_lookup, by = "variant_key") %>%
    mutate(causal = replace_na(causal, FALSE))
}

# one independent +/- flank_bp window per causal variant
build_causal_promoterAI_facets <- function(promoter_data, cs_palette_vec, flank_bp = 20,
                                           freq_lookup = NULL) {
  plot_data <- promoter_data %>%
    mutate(across(c(chrom, ref, alt, cs_id), as.character), pos = as.numeric(pos),
           causal = coalesce(as.logical(causal), FALSE),
           variant_label = paste0(chrom, ":", format(pos, scientific = FALSE, trim = TRUE),
                                   " ", ref, ">", alt))
  causal_variants <- plot_data %>% filter(causal) %>%
    distinct(chrom, pos, ref, alt, cs_id, variant_label) %>% arrange(chrom, pos, ref, alt) %>%
    mutate(finemap_color = unname(cs_palette_vec[cs_id]),
           variant_key   = normalize_variant_id(make_variant_id(chrom, pos, ref, alt)))
  if (!is.null(freq_lookup))
    causal_variants <- causal_variants %>%
      left_join(freq_lookup %>% select(variant_key, freq_label, afc_label), by = "variant_key") %>%
      mutate(strip_label = if_else(is.na(freq_label), variant_label,
                                   paste0(variant_label, "\n", freq_label, "\n", afc_label)))
  else
    causal_variants <- causal_variants %>% mutate(strip_label = variant_label)
  facet_information <- causal_variants %>%
    transmute(chrom,
              facet_id = paste(chrom, format(pos, scientific = FALSE, trim = TRUE), ref, alt, sep = "_"),
              focal_position = pos, focal_ref = ref, focal_alt = alt,
              window_start = pos - flank_bp, window_end = pos + flank_bp,
              cs_id, finemap_color, facet_strip_label = strip_label) %>%
    arrange(chrom, focal_position) %>%
    mutate(facet_strip_label = factor(facet_strip_label, levels = unique(facet_strip_label)))
  facet_plot_data <- plot_data %>%
    inner_join(facet_information, by = "chrom", relationship = "many-to-many") %>%
    filter(pos >= window_start, pos <= window_end) %>%
    mutate(relative_position = pos - focal_position,
           is_facet_causal = pos == focal_position & ref == focal_ref & alt == focal_alt,
           facet_strip_label = factor(facet_strip_label, levels = levels(facet_information$facet_strip_label)))
  attr(facet_plot_data, "facet_information") <- facet_information
  facet_plot_data
}

plot_causal_promoterAI_facets <- function(promoter_data, cs_palette_vec, flank_bp = 20,
    effect_threshold = 0.4, causal_point_size = 4, facet_diamond_size = 3,
    facet_diamond_stroke = 0.8, motif_annotations = NULL, freq_lookup = NULL,
    nucleotide_cols = nucleotide_colors) {

  fpd        <- build_causal_promoterAI_facets(promoter_data, cs_palette_vec, flank_bp, freq_lookup)
  facet_info <- attr(fpd, "facet_information")
  causal_pts <- fpd %>% filter(is_facet_causal) %>%
    distinct(facet_id, facet_strip_label, chrom, pos, ref, alt, promoterAI, relative_position, .keep_all = TRUE)
  facet_symbols <- facet_info %>% transmute(facet_id, facet_strip_label, relative_position = 0, finemap_color)

  y_limit  <- max(1, ceiling(max(abs(fpd$promoterAI), na.rm = TRUE) * 10) / 10)
  y_breaks <- sort(unique(c(-y_limit, -0.75, -0.5, 0, 0.5, 0.75, y_limit)))
  y_breaks <- y_breaks[y_breaks >= -y_limit & y_breaks <= y_limit]

  motif_col    <- "#3E86C0"
  motif_layers <- if (!is.null(motif_annotations)) {
    ann <- motif_annotations %>%
      left_join(facet_info %>% select(facet_id, facet_strip_label), by = "facet_id") %>%
      mutate(facet_strip_label = factor(facet_strip_label, levels = levels(fpd$facet_strip_label)))
    list(
      geom_rect(data = ann, aes(xmin = xmin, xmax = xmax), ymin = -Inf, ymax = y_limit,
                fill = motif_col, alpha = 0.22, inherit.aes = FALSE),
      geom_text(data = ann, aes(x = (xmin + xmax) / 2, label = paste0(tf, "\nmotif")),
                y = y_limit * 0.82, inherit.aes = FALSE, fontface = "bold.italic",
                size = 3.2, colour = "grey15", lineheight = 0.85))
  } else NULL

  ggplot(fpd, aes(relative_position, promoterAI)) +
    motif_layers +
    annotate("segment", x = -flank_bp * 0.88, xend = flank_bp * 0.88,
             y = y_limit * 1.12, yend = y_limit * 1.12,
             colour = "grey72", linewidth = 0.6, lineend = "round") +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.45) +
    geom_hline(yintercept = c(-effect_threshold, effect_threshold),
               linetype = "dotted", color = "grey68", linewidth = 0.4) +
    geom_segment(aes(xend = relative_position, y = 0, yend = promoterAI),
                 color = "grey72", linewidth = 0.35, alpha = 0.75) +
    geom_point(data = filter(fpd, !is_facet_causal, abs(promoterAI) < effect_threshold),
               shape = 16, size = 1.5, color = "grey76", alpha = 0.55, show.legend = FALSE) +
    geom_point(data = filter(fpd, !is_facet_causal, abs(promoterAI) >= effect_threshold),
               aes(color = alt), shape = 16, size = 2.3, alpha = 0.9) +
    geom_point(data = causal_pts, aes(fill = alt), shape = 21, size = causal_point_size,
               stroke = 0.7, color = "black", show.legend = FALSE) +
    geom_point(data = facet_symbols, aes(relative_position, y = Inf, fill = I(finemap_color)),
               inherit.aes = FALSE, shape = 23, size = facet_diamond_size,
               stroke = facet_diamond_stroke, color = "black", show.legend = FALSE) +
    facet_grid(cols = vars(facet_strip_label), scales = "free_x", space = "free_x") +
    scale_color_manual(values = nucleotide_cols, breaks = c("A", "C", "T", "G"),
                       drop = FALSE, name = "Alternate allele") +
    scale_fill_manual(values = nucleotide_cols, breaks = c("A", "C", "T", "G"), guide = "none") +
    scale_x_continuous(breaks = \(x) pretty(x, 5), expand = expansion(mult = c(0.025, 0.025))) +
    scale_y_continuous(trans = signed_square_trans(), limits = c(-y_limit, y_limit),
                       oob = scales::oob_keep, breaks = y_breaks,
                       labels = label_number(accuracy = 0.01, trim = TRUE),
                       expand = expansion(mult = c(0.05, 0.15))) +
    coord_cartesian(clip = "off") +
    labs(x = "Position relative to causal variant (bp)", y = "PromoterAI effect") +
    guides(color = guide_legend(title.position = "top", ncol = 2, byrow = TRUE,
                                override.aes = list(size = 2.5, alpha = 1))) +
    theme(strip.text.x = element_text(size = 9, face = "bold", margin = margin(2, 2, 15, 2)),
          panel.spacing.x = unit(1.2, "lines"),
          axis.text = element_text(size = 9), axis.title = element_text(size = 11),
          legend.position = "inside", legend.position.inside = c(0.02, 0.04),
          legend.justification = c(0, 0), legend.title = element_text(size = 7),
          legend.text = element_text(size = 7), legend.key.size = unit(3.5, "mm"),
          legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
          plot.margin = margin(8, 15, 5, 5))
}

.acvrl1_qtl <- function(eqtl)
  eqtl %>% filter(gene_name == "ACVRL1") %>%
    mutate(causal = pip > 0.9, cs_id = as.character(cs_id))

build_panelD <- function(eqtl) {
  qtl <- .acvrl1_qtl(eqtl)
  plot_finemapping_with_inset(
    qtl %>% mutate(type = "Expression"),
    pip_thresh = 0.9, zoom_center = ACVRL1_FOCAL, zoom_window_bp = 500,
    cs_colors = cs_palette(qtl$cs_id),
    panel_title = "ACVRL1 eQTL fine-mapping", ccre = fread(input_path("acvrl1_ccre")))
}

build_panelE <- function(eqtl, variant_freq) {
  qtl <- .acvrl1_qtl(eqtl)
  promoterAI <- prepare_promoterAI_locus(
    fread(input_path("acvrl1_promoterai")), qtl, gene_symbol = "ACVRL1",
    transcript = "ENST00000419526.6", focal_position = ACVRL1_FOCAL, window_bp = 100)
  # focal green-diamond variant (chr12:51912316 C>T) sits in an SP1 GC-box
  plot_causal_promoterAI_facets(
    promoterAI, cs_palette_vec = cs_palette(qtl$cs_id), flank_bp = 20, effect_threshold = 0.5,
    freq_lookup = variant_freq,
    motif_annotations = tibble(facet_id = "chr12_51912316_C_T", xmin = -3.5, xmax = 5.5, tf = "SP1"))
}

# =============================================================================
# Panels F/G (SPR worked example). Borzoi ref/alt bigWigs come from
# input_path("spr_bigwig") (a directory) and the attribution table from
# input_path("spr_attr"). rtracklayer is used namespaced (never attached).
# Objects are UNTAGGED; tags are applied in fig5/assemble.R.
# =============================================================================
SPR_VARIANT_POS <- 72887998                                       # SPR variant, chr2 (hg38)

MARK_ORDER <- c("CTCF", "cohesin", "ELF1", "accessibility", "DNase", "ATAC", "EP300",
                "H3K27ac", "H3K4me1", "H3K4me3", "CAGE", "RNA", "H3K36me3", "EZH2",
                "H3K27me3", "H3K9me3")

# per-mark reference-track colours (cool tones, red-free so the red alt line contrasts)
MARK_COLORS <- c(
  CTCF = "#253473", cohesin = "#3CC8C0", ELF1 = "#61619B", accessibility = "#67A850",
  DNase = "#67A850", ATAC = "#41AE76", EP300 = "#8CCCE4", H3K27ac = "#4982B6",
  H3K4me1 = "#4E9678", H3K4me3 = "#303251", CAGE = "#A0B8C3", RNA = "#2D332D",
  H3K36me3 = "#B0D9CB", EZH2 = "#6F8BA0", H3K27me3 = "#3E5B59", H3K9me3 = "#4A5F71")

# broad -> narrow progressive zoom; base-resolution motif logo goes last
ZOOM_GROUPS <- list(
  "repressive / broad"   = list(half = 3000, label = "Repressive histone",
      marks = c("H3K36me3", "EZH2", "H3K27me3", "H3K9me3")),
  "active histone"       = list(half = 1500, label = "Activating histone",
      marks = c("EP300", "H3K27ac", "H3K4me1", "H3K4me3", "RNA")),
  "CTCF / accessibility" = list(half = 500,  label = "CTCF / accessibility",
      marks = c("DNase", "cohesin", "CTCF", "CAGE")))

# ENCODE SCREEN cCREs at the SPR locus (hg38), relative to the variant;
# the variant lies inside EH38E2009446 (pELS, CTCF-bound).
ENCODE_CCRE <- tribble(
  ~name,          ~start, ~end,  ~type,
  "EH38E2009442",  -1432, -1169, "pELS",
  "EH38E2009443",  -1096,  -801, "pELS,CTCF-bound",
  "EH38E2009444",   -789,  -492, "PLS,CTCF-bound",
  "EH38E2009445",   -382,  -171, "pELS,CTCF-bound",
  "EH38E2009446",   -163,    70, "pELS,CTCF-bound",
  "EH38E2009447",    538,   867, "pELS,CTCF-bound",
  "EH38E2009448",    905,  1206, "pELS,CTCF-bound")
ENCODE_ELEMENT <- with(dplyr::filter(ENCODE_CCRE, start <= 0, end >= 0), c(start[1], end[1]))

# --- Borzoi bigWig loading (rtracklayer namespaced) ------------------------
.spr_build_metadata <- function(bw_dir, mark_order = MARK_ORDER) {
  tibble(ref_path = list.files(bw_dir, pattern = "\\.ref\\.bw$",
                               recursive = TRUE, full.names = TRUE)) |>
    mutate(alt_path = str_replace(ref_path, "\\.ref\\.bw$", ".alt.bw"),
           mark     = str_remove(basename(ref_path), "\\.ref\\.bw$"),
           mark     = factor(mark, levels = mark_order)) |>
    arrange(mark) |> select(mark, ref_path, alt_path)
}
.spr_read_bw <- function(path, label)
  as.data.frame(rtracklayer::import(path)) |>
    transmute(score, type = factor(label, levels = c("reference", "alternate")), start)
.spr_load_track <- function(mark, meta, variant_pos = SPR_VARIANT_POS) {
  row <- filter(meta, mark == !!mark)
  bind_rows(.spr_read_bw(row$ref_path[1], "reference"),
            .spr_read_bw(row$alt_path[1], "alternate")) |>
    mutate(mark = factor(as.character(mark), levels = levels(meta$mark)),
           relative_pos = start - variant_pos)
}
.spr_load_all <- function(meta, variant_pos = SPR_VARIANT_POS)
  map_dfr(as.character(meta$mark), .spr_load_track, meta = meta, variant_pos = variant_pos)

# % change at the regulatory element (each mark's ref peak +/- half)
pct_at_element <- function(df, search = c(-750, 750), half = 96) {
  d <- filter(df, relative_pos >= search[1], relative_pos <= search[2])
  peaks <- d |> filter(type == "reference") |>
    group_by(mark) |> slice_max(score, n = 1, with_ties = FALSE) |>
    transmute(mark, peak_pos = relative_pos)
  d |> inner_join(peaks, by = "mark") |>
    filter(relative_pos >= peak_pos - half, relative_pos <= peak_pos + half) |>
    group_by(mark, type) |> summarise(sig = mean(score), .groups = "drop") |>
    pivot_wider(names_from = type, values_from = sig) |>
    mutate(pct = 100 * (`alternate` - reference) / reference) |> arrange(mark)
}
.spr_labeller <- function(df, search = c(-750, 750), half = 96) {
  pc <- pct_at_element(df, search, half)
  as_labeller(setNames(sprintf("%s\n(%+.0f%%)", pc$mark, pc$pct), as.character(pc$mark)))
}

# grouped browser (broad -> narrow), cool per-mark fill + red alt line
plot_browser_grouped <- function(df, groups = ZOOM_GROUPS, mark_colors = MARK_COLORS,
                                 alt_col = "#d1495b", fill_alpha = 0.25, ref_lwd = 0.7,
                                 min_pct = 10, drop = c("RNA"),
                                 element = NULL, hl_fill = "#F4A261", hl_alpha = 0.22,
                                 bottom = NULL, bottom_height = 2.5) {
  pc <- pct_at_element(df)
  keep_marks <- setdiff(as.character(pc$mark[abs(pc$pct) >= min_pct]), drop)
  groups <- purrr::keep(map(groups, \(g) { g$marks <- intersect(g$marks, keep_marks); g }),
                        \(g) length(g$marks) > 0)
  df <- filter(df, mark %in% keep_marks)
  if (is.null(element)) element <- ENCODE_ELEMENT

  grp_panel <- function(marks, half, x_title = FALSE) {
    d <- df |> filter(mark %in% marks) |>
         mutate(mark = factor(as.character(mark), levels = marks)) |>
         filter(relative_pos >= -half, relative_pos <= half)
    fmt <- \(v) ifelse(v < 3, sprintf("%.1f", v), sprintf("%.0f", v))
    rng <- d |> group_by(mark) |>
           summarise(lab = sprintf("[0-%s]", fmt(max(score))), .groups = "drop")
    ggplot(d, aes(relative_pos, score)) +
      annotate("rect", xmin = element[1], xmax = element[2], ymin = -Inf, ymax = Inf,
               fill = hl_fill, alpha = hl_alpha) +
      geom_area(data = \(x) filter(x, type == "reference"), aes(fill = mark), alpha = fill_alpha) +
      geom_line(data = \(x) filter(x, type == "reference"), aes(colour = mark), linewidth = ref_lwd) +
      geom_line(data = \(x) filter(x, type == "alternate"), colour = alt_col, linewidth = .6) +
      geom_vline(xintercept = 0, colour = "grey30", linewidth = .4) +
      geom_text(data = rng, aes(x = -half, y = Inf, label = lab), inherit.aes = FALSE,
                hjust = -0.08, vjust = 1.35, size = 2.3, colour = "grey45") +
      scale_fill_manual(values = mark_colors, guide = "none") +
      scale_colour_manual(values = mark_colors, guide = "none") +
      scale_x_continuous(breaks = c(-half, -half / 2, 0, half / 2, half)) +
      facet_grid(rows = vars(mark), scales = "free_y", switch = "y",
                 labeller = .spr_labeller(df |> filter(mark %in% marks))) +
      labs(x = if (x_title) "position relative to variant (bp)" else NULL, y = NULL) +
      theme_classic(base_size = 10) +
      theme(strip.placement = "outside", strip.background = element_blank(),
            strip.text.y.left = element_text(angle = 0, size = 9.5), panel.spacing = unit(2, "pt"),
            axis.text.y = element_blank(), axis.ticks.y = element_blank(),
            axis.line.y = element_blank(), legend.position = "none")
  }
  panels <- map(groups, \(g) grp_panel(g$marks, g$half))
  rel_h  <- map_dbl(groups, \(g) length(g$marks) + 1)
  if (!is.null(bottom)) { panels <- c(panels, list(bottom)); rel_h <- c(rel_h, bottom_height) }
  body <- plot_grid(plotlist = panels, ncol = 1, align = "v", axis = "lr", rel_heights = rel_h)

  ctr_top <- cumsum(rel_h)[seq_along(groups)] - rel_h[seq_along(groups)] / 2
  yf      <- 1 - ctr_top / sum(rel_h)
  right   <- ggdraw()
  for (i in seq_along(groups))
    right <- right + draw_label(map_chr(groups, "label")[i], x = 0.55, y = yf[i],
                                angle = 270, size = 8.5, fontface = "bold", colour = "grey30")
  plot_grid(body, right, ncol = 2, rel_widths = c(1, 0.05))
}

# base-resolution attribution motif logo
motif_logo <- function(attr_file, window = c(-8, 12), target = "CTCF", hl = "#d1495b",
                       hl_extend_right = 0) {
  a     <- read_tsv(attr_file, show_col_types = FALSE) |> filter(rel >= window[1], rel <= window[2])
  bases <- c("A", "C", "G", "T")
  mat   <- matrix(0, 4, nrow(a), dimnames = list(bases, as.character(a$rel)))
  for (i in seq_len(nrow(a))) mat[a$base[i], i] <- a$attr[i] / 1000
  var_cols <- which(a$deleted)
  p <- ggseqlogo(mat, method = "custom")
  if (length(var_cols))
    p <- p + annotate("rect", xmin = min(var_cols) - 0.5, xmax = max(var_cols) + 0.5 + hl_extend_right,
                      ymin = -Inf, ymax = Inf, fill = hl, alpha = 0.18)
  p +
    annotate("text", x = 0.5, y = Inf, label = sprintf("[%.0f, %.0f]", min(mat), max(mat)),
             hjust = 0, vjust = 1.35, size = 2.3, colour = "grey45") +
    scale_x_continuous(breaks = seq(1, nrow(a), 4), labels = a$rel[seq(1, nrow(a), 4)]) +
    labs(x = "position relative to variant (bp)", y = sprintf("%s\nattribution", target)) +
    theme_classic(base_size = 10) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
          axis.title.y = element_text(size = 8, colour = "grey20", angle = 0, vjust = 0.5))
}

# cowplot browser: grouped tracks + CTCF attribution logo at the bottom
make_spr_browser <- function(bw_dir = input_path("spr_bigwig"), attr_file = input_path("spr_attr")) {
  meta <- .spr_build_metadata(bw_dir)
  # ref allele GCGCC spans 5 bp but only 4 are flagged deleted; extend one base right
  logo <- motif_logo(attr_file, target = "CTCF", hl_extend_right = 1)
  plot_browser_grouped(.spr_load_all(meta), bottom = logo, bottom_height = 3.6)
}

# --- panel builders (untagged) ---------------------------------------------
build_panelF <- function(eqtl) {
  spr_region <- eqtl %>%
    filter(gene_name == "SPR", between(position / 1e6, 72.881, 72.8925)) %>%
    mutate(causal = pip > 0.9, cs_id = as.character(cs_id), type = "Expression")
  plot_finemapping_with_inset(
    spr_region, pip_thresh = 0.9,
    zoom_center = 72887978, zoom_window_bp = 200, x_limits = c(72.8817, 72.8934),
    cs_colors = cs_palette(spr_region$cs_id),
    panel_title = "SPR eQTL fine-mapping", ccre = fread(input_path("spr_ccre")))
}

build_panelG <- function(eqtl, variant_freq) {
  spr_qtl <- eqtl %>% filter(gene_name == "SPR") %>%
    mutate(causal = pip > 0.9, cs_id = as.character(cs_id), is_focal = position == SPR_VARIANT_POS)
  pal <- cs_palette(spr_qtl$cs_id)
  focal <- spr_qtl %>% filter(is_focal) %>%
    mutate(variant_key = normalize_variant_id(variant)) %>%
    left_join(variant_freq %>% select(variant_key, freq_label, afc_label), by = "variant_key") %>%
    dplyr::slice(1)
  focal_color <- unname(pal[as.character(focal$cs_id)])
  focal_vid   <- sprintf("chr2:%s %s>%s",
                         format(focal$position, big.mark = ",", trim = TRUE), focal$ref, focal$alt)
  # header carries the focal-variant title + credible-set diamond (the panel tag
  # itself is applied in fig5/assemble.R, not here).
  g_header <- ggplot() +
    annotate("text", x = 0.5, y = 0.94, hjust = 0.5, vjust = 1, lineheight = 0.9,
             label = paste(focal_vid, focal$freq_label, focal$afc_label, sep = "\n"),
             size = 3.16, fontface = "bold", family = FONT, colour = "grey15") +
    annotate("point", x = 0.5, y = 0.16, shape = 23, size = 3.6,
             fill = focal_color, colour = "black", stroke = 0.8) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_void(base_family = FONT)
  wrap_elements(full = cowplot::plot_grid(g_header, make_spr_browser(),
                                          ncol = 1, rel_heights = c(0.14, 1)))
}
