# =============================================================================
# Figure 1 — shared setup (sourced by every prep / panel / assembly script)
# -----------------------------------------------------------------------------
# Run all scripts from the PROJECT ROOT (the folder containing input_data/ and
# fig1/). Each script starts with:  source("fig1/00_setup.R")
# Re-sourcing is cheap and idempotent.
#
# The analysis / plotting functions below were folded in from
# QTLAnalysisUtils/QTLGeneAnalysis (LoadExternalQTLs.R, QTLGeneAnalysisFunctions.R,
# PlotFunctions.R) so this figure is self-contained. Bugs present in the originals
# are fixed here (see docs/superpowers/specs/2026-07-23-fig1-cleanup-design.md).
# =============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(rtracklayer)
  library(scales)
  library(ggrepel)
  library(patchwork)
  library(readxl)
  library(showtext)
  library(sysfonts)
  library(paletteer)
  library(rcartocolor)
  library(colorspace)
})

# --- Paths (override via env vars if needed) --------------------------------
INPUT    <- Sys.getenv("FIG1_INPUT",  "input_data")   # raw inputs
DATA_DIR <- Sys.getenv("FIG1_DATA",   "fig1/data")    # per-panel data (01_prep)
OBJ_DIR  <- Sys.getenv("FIG1_OBJ",    "fig1/objects") # per-panel ggplot objects
OUTDIR   <- Sys.getenv("FIG1_OUTDIR", ".")            # final Figure1.pdf/png
for (d in c(DATA_DIR, OBJ_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# Convenience: build an input path (optionally under a subdirectory).
in_path <- function(...) file.path(INPUT, ...)

# --- Font: Inter via Google if reachable, else default sans -----------------
FONT <- "sans"
if (!"inter" %in% sysfonts::font_families()) {
  if (isTRUE(tryCatch({ font_add_google("Inter", "inter"); TRUE }, error = function(e) FALSE)))
    FONT <- "inter"
} else FONT <- "inter"

DPI <- 300                     # showtext DPI; must match every ggsave() DPI
showtext_opts(dpi = DPI)
showtext_auto()

# --- Design system: uniform text sizes --------------------------------------
# element_text() sizes are in pt; geom_text()/annotate() sizes are in mm.
SZ_AXIS_TITLE <- 25; SZ_AXIS_TEXT <- 20; SZ_STRIP <- 20
SZ_LEG_TEXT   <- 15; SZ_LEG_TITLE <- 17; SZ_TAG   <- 30

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

# --- Ancestry colour palette (shared across panels) -------------------------
AncestryPalette <- carto_pal(8, "Temps")
AncestryColorPalette <- data.frame(
  Group = c("ADMIX", "EUR", "AMR", "AFR", "SAS", "EAS", "MID", "AoU"),
  Color = AncestryPalette
)

# =============================================================================
# Folded-in functions (from QTLAnalysisUtils, bug-fixed)
# =============================================================================

# --- External-QTL loaders (LoadExternalQTLs.R) ------------------------------
LoadIntervaleQTL <- function(path) {
  read_xlsx(path, skip = 3) %>%
    mutate(Significant = cis_eQTL_significant == 1) %>%
    select(feature_id, Significant) %>%
    dplyr::rename("gene_id" = "feature_id") %>%
    mutate(study = "INTERVAL")
}

LoadIntervalsQTL <- function(path) {
  read_xlsx(path, skip = 3) %>%
    select(gene_id) %>%
    mutate(gene_id = str_extract(gene_id, "ENSG[0-9]+")) %>%
    mutate(Significant = TRUE, study = "INTERVAL")
}

LoadGTExeQTL <- function(path) {
  fread(path) %>%
    mutate(Significant = qval < .05, gene_id = str_remove(gene_id, "\\..*")) %>%
    select(gene_id, Significant) %>%
    mutate(study = "GTEx")
}

LoadGTExsQTL <- function(path) {
  fread(path) %>%
    mutate(Significant = qval < .05, gene_id = str_remove(gene_id, "\\..*")) %>%
    select(gene_id, Significant) %>%
    mutate(study = "GTEx")
}

LoadTopmedeQTL <- function(path) {
  fread(path) %>%
    mutate(Significant = qvalue < .05, gene_id = str_remove(phenotype_id, "\\..*")) %>%
    select(gene_id, Significant) %>%
    mutate(study = "Topmed")
}

# UKBB pQTL reference. Currently deferred (no input file); kept for the
# novel-pGene comparison once the reference is supplied. See NOTES.md.
LoadUKBBpQTL <- function(path) {
  read_xlsx(path, sheet = 10, skip = 4) %>%
    filter(`cis/trans` == "cis") %>%
    dplyr::rename("gene_id" = `Ensembl gene ID`)
}

LoadExternaleQTLs <- function(GTExPath, IntervalPath, TopmedPath) {
  bind_rows(
    LoadGTExeQTL(GTExPath),
    LoadIntervaleQTL(IntervalPath),
    LoadTopmedeQTL(TopmedPath)
  )
}

LoadExternalsQTLs <- function(GTExPath, IntervalPath) {
  # FIX: original bound the undefined GTEx_eQTL_df; bind the sQTL frame.
  bind_rows(
    LoadGTExsQTL(GTExPath),
    LoadIntervalsQTL(IntervalPath)
  )
}

# --- QTL summarisation / transforms (QTLGeneAnalysisFunctions.R) ------------
ExtractNovelQTLs <- function(input_QTLs, QTL_list) {
  input_QTLs %>%
    filter(!phenotype_id %in% QTL_list) %>%
    dplyr::count(gene_type) %>%
    filter(!str_detect(gene_type, "pseudo|IG|TEC|TR")) %>%
    mutate(group = case_when(
      gene_type == "lncRNA" ~ "lncRNA",
      gene_type == "protein_coding" ~ "Protein-Coding",
      TRUE ~ "non-coding RNA"
    )) %>%
    group_by(group) %>%
    summarize(count = sum(n))
}

SummarizeSigQTLs <- function(QTL_data, AncestryCounts) {
  QTL_data %>%
    mutate(group = toupper(group)) %>%
    group_by(group, Assay) %>%
    mutate(sig = FDR < .05) %>%
    summarize(prop_sig = mean(sig), num_sig = dplyr::n()) %>%
    # FIX: join the AncestryCounts argument, not a global.
    left_join(AncestryCounts, by = c("group" = "ancestry_pred_other")) %>%
    mutate(group = as.factor(group)) %>%
    mutate(shape = ifelse(group == "COMB", 23, 21)) %>%
    mutate(label = case_when(group == "COMB" ~ paste0("COMBINED", "\n(", num_sig, ")"), TRUE ~ group)) %>%
    mutate(group = case_when(group == "COMB" ~ "AoU", TRUE ~ group))
}

ComputeCredibleSetCDF <- function(cs_count_data) {
  cs_count_data %>%
    group_by(group, Assay, molecular_trait_id) %>%
    summarize(number_cs = dplyr::n()) %>%
    group_by(group, number_cs) %>%
    summarise(n = dplyr::n(), .groups = "drop") %>%
    complete(group, number_cs = 0:max(number_cs), fill = list(n = 0)) %>%
    group_by(group) %>%
    arrange(number_cs) %>%
    mutate(prop_ge = rev(cumsum(rev(n))) / sum(n)) %>%   # survival curve
    ungroup()
}

# --- Fine-mapping annotation enrichment (PlotFunctions.R) -------------------
ComputeFinemappingAnnotationEnrichment <- function(AnnotatedFinemapping) {
  VariantAnnotationsLong <- AnnotatedFinemapping %>%
    filter(pip > .9 | pip < .01) %>%
    select(variant, where(is.logical), pip) %>%
    group_by(variant) %>%
    filter(pip == max(pip)) %>%
    mutate(causal = pip > .9) %>%
    select(-pip) %>%
    ungroup() %>%
    distinct() %>%
    pivot_longer(!c(variant, causal)) %>%
    na.omit()
  VariantAnnotationsLong %>%
    group_by(name) %>%
    filter(
      dplyr::n_distinct(causal) == 2,
      dplyr::n_distinct(value) == 2
    ) %>%
    do({
      ft <- fisher.test(table(.$causal, .$value))
      tibble(
        odds_ratio = unname(ft$estimate),
        conf.low   = ft$conf.int[1],
        conf.high  = ft$conf.int[2],
        p.value    = ft$p.value
      )
    }) %>%
    ungroup()
}

# =============================================================================
# Folded-in plot builders (from PlotFunctions.R)
# =============================================================================

PlotAncestryCounts <- function(AncestryCounts) {
  AncestryCounts %>%
    left_join(AncestryColorPalette, by = c("ancestry_pred_other" = "Group")) %>%
    ggplot(aes(x = reorder(toupper(ancestry_pred_other), -n), y = n, fill = Color)) +
    geom_col(show.legend = FALSE, lwd = 1) +
    ylab("Number of participants") +
    xlab("Ancestry Group") +
    theme(axis.title.y = element_text(margin = margin(r = 20.5)),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank()) +
    scale_fill_identity() +
    scale_color_identity()
}

PlotNovelQTLCount <- function(plot_df) {
  pd <- position_dodge(width = 0.8)
  ggplot(plot_df, aes(x = group, y = count + 1, fill = type, group = type)) +
    geom_col(position = pd, width = 0.75, na.rm = TRUE) +
    geom_text(data = dplyr::filter(plot_df, !is.na(count)),
              aes(label = count), position = pd, vjust = -0.2, size = 5) +
    scale_y_log10() +
    xlab("Gene Type") +
    ylab("Novel QTL genes") +
    paletteer::scale_fill_paletteer_d("ltc::minou") +
    paletteer::scale_color_paletteer_d("ltc::minou") +
    theme(legend.position = c(0.93, 0.85),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank())
}

PlotPropQTLSig <- function(plot_df) {
  plot_df %>%
    ggplot(aes(x = n / 1000, y = prop_sig)) +
    geom_point(aes(shape = group, fill = group),
               size = 12, color = "black", stroke = .5) +
    scale_shape_manual(values = c(
      AFR = 21, AMR = 21, AoU = 23,  # 23 = diamond
      EAS = 21, EUR = 21, MID = 21, SAS = 21
    )) +
    facet_wrap(~Assay, scales = "free", nrow = 2) +
    scale_fill_manual(values = c(AFR = "#C8D892",
                                 AMR = "#7FC385",
                                 AoU = "#CF597E",
                                 EAS = "#ECA676",
                                 EUR = "#30AC86",
                                 MID = "#E47D72",
                                 SAS = "#EBCE8D")) +
    guides(fill = guide_legend(override.aes = list(
      size = 8, color = "black", stroke = 1
    ))) +
    ylab("Proportion of genes \nwith significant QTL") +
    xlab("Sample Size (Thousands)") +
    scale_x_log10(breaks = c(0.3, 1, 3, 10),
                  labels = c("0.3k", "1k", "3k", "10k"),
                  expand = expansion(mult = c(0.1, 0.1))) +
    scale_y_continuous(expand = expansion(mult = c(0.15, 0.15)))
}

PlotCSSizes <- function(plot_data, AncestryColorPalette) {
  cleaned_data <- plot_data %>%
    mutate(group = case_when(group == "COMB" ~ "AoU", TRUE ~ group)) %>%
    left_join(AncestryColorPalette, by = c("group" = "Group")) %>%
    mutate(
      Color_fill = lighten(Color, 0),
      Color_line = darken(Color, 0.4)   # darker outline
    )
  cleaned_data %>%
    ggplot(aes(x = fct_reorder(group, cs_size, .fun = mean), y = cs_size,
               fill = Color_fill, color = Color_line)) +
    geom_boxplot(show.legend = FALSE, lwd = 0.4) +
    scale_fill_identity() +
    scale_color_identity() +
    xlab("Ancestry Group") +
    scale_y_log10() +
    ylab("Credible set size")
}

# FIX: takes the merged enrichment data and the functional-group table as
# arguments instead of reading module-level globals.
PlotFinemappingEnrichment <- function(MergedEnrichmentData, variant_functional_groups) {
  annotation_labels <- c(
    splice_donor_variant    = "Splice donor",
    splice_acceptor_variant = "Splice acceptor",
    stop_gained             = "Stop gained",
    missense_variant        = "Missense",
    inframe_insertion       = "Inframe insertion",
    inframe_deletion        = "Inframe deletion",
    PLS                     = "Promoter",
    dELS                    = "Distal enhancer",
    `5_prime_utr_variant`   = "5' UTR",
    `3_prime_utr_variant`   = "3' UTR",
    pELS                    = "Proximal enhancer",
    splice_region_variant   = "Splice region"
  )
  plot_data <- MergedEnrichmentData %>%
    mutate(padj = p.adjust(p.value, method = "fdr")) %>%
    filter(name != "CTCF.only") %>%
    filter(padj < .05) %>%
    left_join(variant_functional_groups, by = "name") %>%
    filter(name %in% names(annotation_labels)) %>%
    mutate(annotation_clean = annotation_labels[name]) %>%
    group_by(functional_group) %>%
    mutate(row = as.numeric(factor(annotation_clean))) %>%
    ungroup()
  shade_df <- plot_data %>%
    distinct(functional_group, annotation_clean, row) %>%
    filter(row %% 2 == 0)
  xmax_rect <- max(plot_data$conf.high, na.rm = TRUE)
  xmin_rect <- min(1, min(plot_data$conf.low, na.rm = TRUE) / 1.5)
  xmax_rect <- xmax_rect * 1.5

  plot_data %>%
    mutate(annotation_clean = case_when(
      str_detect(annotation_clean, "UTR") ~ annotation_clean,
      TRUE ~ str_replace(annotation_clean, " ", "\n")
    )) %>%
    ggplot(aes(y = odds_ratio, x = reorder(annotation_clean, odds_ratio), fill = type)) +
    geom_rect(data = shade_df,
              aes(xmin = row - .5, xmax = row + .5,
                  ymin = xmin_rect - .1, ymax = xmax_rect + .1),
              inherit.aes = FALSE, fill = "grey90") +
    geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                    color = "black", show.legend = FALSE,
                    position = position_dodge(width = 1)) +
    geom_point(size = 7, shape = 21, position = position_dodge(width = 1)) +
    facet_wrap(~functional_group, scales = "free_x", nrow = 1, strip.position = "top") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_y_continuous(trans = "log10", expand = expansion(mult = c(0.015, 0.015))) +
    xlab("") +
    ylab("Fold Enrichment") +
    theme(axis.text.x = element_text(size = 12),
          legend.position = "bottom") +
    paletteer::scale_fill_paletteer_d("ltc::minou")
}

LoadAncestryCounts <- function(path, sample_list = NULL) {
  ancestry_df <- fread(path) %>% select(research_id, ancestry_pred_other)
  if (!is.null(sample_list)) {
    ancestry_df <- ancestry_df %>% filter(research_id %in% sample_list)
  }
  ancestry_df %>%
    dplyr::count(ancestry_pred_other) %>%
    mutate(ancestry_pred_other = case_when(ancestry_pred_other == "oth" ~ "ADMIX",
                                           TRUE ~ ancestry_pred_other)) %>%
    mutate(ancestry_pred_other = toupper(ancestry_pred_other))
}
