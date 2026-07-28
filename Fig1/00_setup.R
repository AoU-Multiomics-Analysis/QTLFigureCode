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
# INPUT is the staging root; it defaults to the current working directory, so
# gsutil-staged files land under the cwd (in their manifest subpaths, e.g.
# permutations/eqtl/…). Set FIG1_INPUT to keep inputs in a separate folder.
INPUT    <- Sys.getenv("FIG1_INPUT",    ".")            # staging root (cwd)
OBJ_DIR  <- Sys.getenv("FIG1_OBJ",      "fig1/objects") # per-panel ggplot objects
OUTDIR   <- Sys.getenv("FIG1_OUTDIR", "fig1")         # final Figure1.pdf/png
MANIFEST <- Sys.getenv("FIG1_MANIFEST", "fig1/inputs.tsv")
dir.create(OBJ_DIR, showWarnings = FALSE, recursive = TRUE)

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
SZ_LEG_TEXT   <- 15; SZ_LEG_TITLE <- 17; SZ_TAG   <- 31

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

# --- Assay colour palette (shared by the assay-coloured panels) -------------
# Keyed by name, not by factor-level order. The novel-gene panel lists its bars
# eGene / sGene / pGene while every other assay-coloured panel is ordered
# expression / protein / splicing, so taking colours positionally off the
# palette gave the splicing bar the protein colour and vice versa. Every alias
# for one assay maps to one colour here, so panel order can never move a colour
# again.
AssayPalette <- as.character(paletteer::paletteer_d("ltc::minou"))[1:3]
AssayColors <- c(
  Expression = AssayPalette[1], eQTL = AssayPalette[1], eGene = AssayPalette[1],
  Protein    = AssayPalette[2], pQTL = AssayPalette[2], pGene = AssayPalette[2],
  Splicing   = AssayPalette[3], sQTL = AssayPalette[3], sGene = AssayPalette[3]
)

# =============================================================================
# Input manifest & lazy staging (shared implementation)
# -----------------------------------------------------------------------------
# fig1/inputs.tsv columns: key, group, Assay, local_subpath, gs_uri
# `key`, `local_subpath` and `gs_uri` are reserved; `group` and `Assay` are
# metadata attached to the rows returned by load_input(). Each script loads only
# the keys it needs (see the panel scripts).
# =============================================================================
dir.create(OUTDIR,   showWarnings = FALSE, recursive = TRUE)

source("common/manifest.R")
source("common/vep.R")

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

# UKBB-PPP cis-pQTL reference for the novel-pGene comparison in Panel D.
# Reads sheet 10 (ST9) of the Nature 2023 supplementary workbook, skipping the
# 4 header rows, keeps the cis associations and exposes the Ensembl gene ID as
# `gene_id`. See NOTES.md.
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
    expand_vep_consequences() %>%
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
    # No colour aesthetic is mapped here, so this panel needs a fill scale only;
    # the colour scale that used to sit alongside it was vestigial.
    scale_fill_manual(values = AssayColors) +
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
    scale_y_log10(labels = scales::label_log()) +
    ylab("Credible set size")
}

# The annotations panel H shows, and how they are labelled. Membership of this
# vector is what puts an annotation on the panel: anything else is filtered out.
#
# VEP spells the UTR terms with a capital UTR. The lowercase form silently
# failed the `name %in% names(...)` filter, dropping both UTR entries from the
# panel even though their enrichments were computed fine (see NOTES.md).
PANEL_ANNOTATION_LABELS <- c(
  splice_donor_variant    = "Splice donor",
  splice_acceptor_variant = "Splice acceptor",
  stop_gained             = "Stop gained",
  missense_variant        = "Missense",
  # The coding slots were inframe insertion and inframe deletion. Frameshift and
  # start lost replaced them: both are significant in all three assays, whereas
  # inframe insertion was not significant for sQTL.
  frameshift_variant      = "Frameshift",
  start_lost              = "Start lost",
  PLS                     = "Promoter",
  dELS                    = "Distal enhancer",
  `5_prime_UTR_variant`   = "5' UTR",
  `3_prime_UTR_variant`   = "3' UTR",
  pELS                    = "Proximal enhancer",
  splice_region_variant   = "Splice region"
)

# NOTE: panel H is now built by PlotGlmEnrichment() below, from the PIP logistic
# regression. This Fisher odds-ratio version is what produced the panel before
# that, and is kept because NOTES.md documents the comparison between the two
# and `fig1/enrichment_glm.R` cross-checks against it. Nothing in the figure
# calls it.
PlotFinemappingEnrichment <- function(MergedEnrichmentData, variant_functional_groups) {
  annotation_labels <- PANEL_ANNOTATION_LABELS
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
  RenderFinemappingEnrichment(plot_data)
}

# Pure rendering seam: takes the prepared per-annotation enrichment table
# (name, type, odds_ratio, conf.low, conf.high, functional_group, row) and
# returns the panel. Split out from PlotFinemappingEnrichment() so the panel can
# be redrawn from an already-computed enrichment table.
RenderFinemappingEnrichment <- function(plot_data) {
  # Fisher's exact test returns an infinite odds ratio, and an infinite upper
  # confidence bound, whenever a cell of the 2x2 table is empty — which happens
  # for the sparsest pQTL annotations. Drawn as-is those points land at +Inf and
  # the panel edge clips them in half. Put them at a finite ceiling above every
  # bounded estimate and draw them as an upward triangle instead, so an
  # unbounded estimate still reads as unbounded.
  bounds <- c(plot_data$odds_ratio, plot_data$conf.high)
  y_cap  <- max(bounds[is.finite(bounds)]) * 1.5
  plot_data <- plot_data %>%
    mutate(
      unbounded   = !is.finite(odds_ratio),
      or_plot     = ifelse(is.finite(odds_ratio), odds_ratio, y_cap),
      high_plot   = ifelse(is.finite(conf.high),  conf.high,  y_cap),
      point_shape = ifelse(unbounded, 24, 21)
    )

  shade_df <- plot_data %>%
    distinct(functional_group, annotation_clean, row) %>%
    filter(row %% 2 == 0)

  # Bands run the full height of the panel (ymax = Inf); the lower edge is
  # finite so it keeps the y = 1 reference line inside the axis range.
  ymin_rect <- min(1, min(plot_data$conf.low, na.rm = TRUE) / 1.5)

  plot_data %>%
    mutate(annotation_clean = case_when(
      str_detect(annotation_clean, "UTR") ~ annotation_clean,
      TRUE ~ str_replace(annotation_clean, " ", "\n")
    )) %>%
    ggplot(aes(y = or_plot, x = reorder(annotation_clean, or_plot), fill = type)) +
    geom_rect(data = shade_df,
              aes(xmin = row - .5, xmax = row + .5,
                  ymin = ymin_rect - .1, ymax = Inf),
              inherit.aes = FALSE, fill = "grey90") +
    geom_pointrange(aes(ymin = conf.low, ymax = high_plot),
                    color = "black", show.legend = FALSE,
                    position = position_dodge(width = 1)) +
    geom_point(aes(shape = point_shape), size = 7,
               position = position_dodge(width = 1)) +
    scale_shape_identity() +
    facet_wrap(~functional_group, scales = "free_x", nrow = 1, strip.position = "top") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    # Headroom at the top so the size-7 markers at the ceiling are not clipped.
    scale_y_continuous(trans = "log10", labels = scales::label_log(),
                       expand = expansion(mult = c(0.015, 0.06))) +
    xlab("") +
    ylab("Fold Enrichment") +
    theme(axis.text.x = element_text(size = 12),
          legend.position = "bottom") +
    scale_fill_manual(values = AssayColors) +
    # shape is an identity aesthetic now, so the fill keys need it spelled out.
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 7)))
}

# --- panel H: annotation enrichment as a PIP logistic regression -------------
# Takes the tidy glm output from FitAnnotationEnrichmentGlm() and the
# functional-group table, and returns the panel. The effect is a log-odds slope,
# so the axis is LINEAR and centred on zero -- there is no odds ratio to log,
# and none of the infinities that a 2x2 table with an empty cell produces.
PlotGlmEnrichment <- function(GlmEnrichmentResults, variant_functional_groups) {
  plot_data <- GlmEnrichmentResults %>%
    filter(name %in% names(PANEL_ANNOTATION_LABELS)) %>%
    left_join(variant_functional_groups, by = "name") %>%
    mutate(annotation_clean = PANEL_ANNOTATION_LABELS[name]) %>%
    group_by(functional_group) %>%
    mutate(row = as.numeric(factor(annotation_clean))) %>%
    ungroup()

  missing_group <- plot_data %>% filter(is.na(functional_group))
  if (nrow(missing_group))
    stop("panel annotations with no functional group: ",
         paste(unique(missing_group$name), collapse = ", "))

  # A separated fit has no finite maximum-likelihood slope, so its interval
  # would swamp the axis. None of the panel annotations separate on the current
  # data; if one ever does, drop it loudly rather than silently rescaling.
  separated <- plot_data %>% filter(separated)
  if (nrow(separated)) {
    message("Panel H: dropping separated fits with no finite slope: ",
            paste(separated$type, separated$name, collapse = ", "))
    plot_data <- plot_data %>% filter(!separated)
  }

  RenderGlmEnrichment(plot_data)
}

# Pure rendering seam, mirroring RenderFinemappingEnrichment().
#
# glm() returns a slope on the log-odds scale; the panel shows it exponentiated,
# on a log10 axis, so the effect reads as an odds ratio between PIP = 1 and
# PIP = 0 and the axis matches the rest of the figure's log panels. exp() is
# monotone, so the confidence bounds transform directly and the ordering of the
# annotations is unaffected.
RenderGlmEnrichment <- function(plot_data) {
  plot_data <- plot_data %>%
    mutate(significant = padj < 0.05,
           or          = exp(estimate),
           or_low      = exp(conf.low),
           or_high     = exp(conf.high))

  # With every annotation significant, the significance key is a legend entry
  # that never varies. It has to be suppressed in guides() as well as on the
  # scale, because an explicit guides() entry overrides the scale's own.
  show_significance <- dplyr::n_distinct(plot_data$significant) > 1

  # A grey rule under each annotation, standing in for the x axis line (the
  # idiom from Fig2 panel B). It spans the width the three dodged assay points
  # occupy, so it reads as "these points belong to this category". `row` runs
  # 1..n within each facet, exactly the x positions reorder() produces, so every
  # category gets one wherever it happens to be ordered.
  # Below everything, including the odds-ratio-of-1 reference line. Taking the
  # minimum over the data alone put the rules on top of that line whenever no
  # annotation was depleted, which is the usual case.
  rule_y <- min(c(plot_data$or_low, 1), na.rm = TRUE) * 0.55
  rule_df <- plot_data %>%
    distinct(functional_group, annotation_clean, row) %>%
    mutate(x = row - 0.4, xend = row + 0.4, y = rule_y, yend = rule_y)

  plot_data %>%
    mutate(annotation_clean = case_when(
      str_detect(annotation_clean, "UTR") ~ annotation_clean,
      TRUE ~ str_replace(annotation_clean, " ", "\n")
    )) %>%
    ggplot(aes(y = or, x = reorder(annotation_clean, or),
               fill = type, alpha = significant)) +
    geom_linerange(aes(ymin = or_low, ymax = or_high),
                   color = "black", show.legend = FALSE,
                   position = position_dodge(width = 1)) +
    geom_point(size = 7, shape = 21, position = position_dodge(width = 1)) +
    geom_segment(data = rule_df, aes(x = x, xend = xend, y = y, yend = yend),
                 inherit.aes = FALSE, color = "grey60", linewidth = 0.7) +
    facet_wrap(~functional_group, scales = "free_x", nrow = 1, strip.position = "top") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    # An odds ratio of 1 is the no-enrichment line, so it stays on the axis even
    # if every interval sits above it.
    expand_limits(y = 1) +
    scale_y_log10(labels = scales::label_log(),
                  expand = expansion(mult = c(0.04, 0.06))) +
    scale_alpha_manual(values = c(`FALSE` = 0.3, `TRUE` = 1), name = NULL,
                       labels = c(`FALSE` = "FDR ≥ 0.05", `TRUE` = "FDR < 0.05")) +
    xlab("") +
    # Two lines: at the figure's axis-title size this runs nearly the full panel
    # height on one.
    ylab("Fine-mapping Enrichment\n(Log Odds)") +
    theme(axis.text.x = element_text(size = 12),
          # The per-category grey rules replace the axis line.
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "bottom") +
    scale_fill_manual(values = AssayColors) +
    # Assay is the primary encoding, so it leads the legend row.
    guides(fill  = guide_legend(order = 1, override.aes = list(alpha = 1, size = 7)),
           alpha = if (show_significance)
                     guide_legend(order = 2, override.aes = list(fill = "grey30", size = 7))
                   else "none")
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
