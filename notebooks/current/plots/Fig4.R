# %% Figure 4 plot theme
#| Construct every plot from the result tables in the current R session.
suppressPackageStartupMessages(library(ggh4x))

fig4_base_theme <- function() {
  theme(
    line = element_line(linewidth = 0.25),
    axis.line = element_line(linewidth = 0.25),
    axis.ticks = element_line(linewidth = 0.25),
    text = element_text(family = "Helvetica", size = 7),
    axis.text = element_text(size = 6, colour = "black"),
    axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6),
    axis.title = element_text(size = 7), strip.text = element_text(size = 7),
    strip.text.x = element_text(size = 7), legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    plot.tag = element_text(size = 10, face = "bold", family = "Helvetica"),
    plot.tag.position = c(0, 1), plot.margin = margin(2, 2, 2, 2, unit = "pt"),
    plot.title = element_blank(), plot.subtitle = element_blank()
  )
}

# %% Trait labels and category colors
#| Preserve the current trait labels and category colors.
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

# --- Short labels for the in-panel trait annotations (Panel A, labelled variant)
# An OVERRIDE map, not a whitelist: only traits whose upstream name is too long,
# too clinical, or an acronym that must not be sentence-cased need an entry.
# Everything else falls through to trait_label()'s tidy-up below. Keys are
# LOWER-CASE and matched case-insensitively — the coloc table has changed the
# case of its trait names before ("Asthma" -> "asthma"), which silently dropped
# every override when the keys were exact. Keep entries to ~15 characters: the
# widest label sets the width of the label column, and that column has to clear
# the longest bar.
TRAIT_SHORTHAND <- c(
  # allergic_respiratory
  "allergic rhinitis" = "Rhinitis", "atopic dermatitis" = "Atopic derm.",
  # anthropometric
  "body mass index" = "BMI",
  # autoimmune
  "inflammatory bowel disease" = "IBD", "rheumatoid arthritis" = "RA",
  "crohn's disease" = "Crohn's", "autoimmune thyroid disease" = "AITD",
  "ulcerative colitis" = "UC", "type 1 diabetes" = "T1D",
  "primary biliary cholangitis" = "PBC", "autoimmune disease" = "Autoimmune",
  "systemic lupus erythematosus" = "SLE", "multiple sclerosis" = "MS",
  # cancer
  "prostate carcinoma" = "Prostate ca.", "breast carcinoma" = "Breast ca.",
  # cardiometabolic
  "type 2 diabetes" = "T2D", "atrial fibrillation" = "AFib",
  "hypertension" = "HTN", "coronary artery disease" = "CAD",
  "abdominal aortic aneurysm" = "AAA", "venous thromboembolism" = "VTE",
  # clinical_chemistry — the assay acronyms would otherwise sentence-case to
  # "Ast" / "Alt" / "Bun".
  "c-reactive protein measurement" = "CRP", "urate measurement" = "Urate",
  "ast" = "AST", "alt" = "ALT", "bun" = "BUN",
  # gastrointestinal
  "cholelithiasis" = "Gallstones", "diverticular disease" = "Diverticular",
  "gastroesophageal reflux disease" = "GERD",
  # hematological — likewise "Mchc" / "Rdw" / "Rbc count".
  "mchc" = "MCHC", "mch" = "MCH", "rdw" = "RDW", "rbc_count" = "RBC count",
  # lipid
  "hdl cholesterol" = "HDL-C", "ldl cholesterol" = "LDL-C",
  "triglycerides" = "TG",
  # metabolic
  "hba1c" = "HbA1c",
  # musculoskeletal
  "hip osteoarthritis" = "Hip OA", "osteoarthritis" = "OA",
  # neurodegenerative
  "alzheimer disease" = "Alzheimer's", "parkinson disease" = "Parkinson's",
  # ophthalmological
  "primary open angle glaucoma" = "POAG",
  "age-related macular degeneration" = "AMD",
  # reproductive_endocrine
  "uterine fibroid" = "Fibroids"
)

# Display label for a raw trait string: the TRAIT_SHORTHAND override if there is
# one (matched case-insensitively), else the upstream name tidied — underscores
# to spaces, and the first letter capitalised only when the name is entirely
# lower-case, so a mixed-case upstream name is left as its authors wrote it.
trait_label <- function(trait) {
  override <- unname(TRAIT_SHORTHAND[tolower(trait)])
  tidied <- gsub("_", " ", trait, fixed = TRUE)
  all_lower <- !is.na(tidied) & tidied == tolower(tidied)
  tidied[all_lower] <- sub("^(.)", "\\U\\1", tidied[all_lower], perl = TRUE)
  ifelse(is.na(override), tidied, override)
}

complete_category_colors <- function(categories) {
  categories <- sort(unique(as.character(categories)))
  missing <- setdiff(categories, names(GROUP_COLORS))
  fallback <- if (length(missing) > 0)
    setNames(scales::hue_pal(l = 55, c = 90)(length(missing)), missing) else character()
  c(GROUP_COLORS, fallback)
}

blend_with_white <- function(hex, amount = 0.65) {
  x <- grDevices::col2rgb(hex) / 255
  y <- x * (1 - amount) + amount
  grDevices::rgb(y[1, ], y[2, ], y[3, ])
}

make_fill <- function(category, intensity, colors) {
  base  <- unname(colors[as.character(category)])
  light <- blend_with_white(base, 0.72)
  vapply(seq_along(base), function(i) {
    palette <- grDevices::colorRampPalette(c(light[[i]], base[[i]]))(100)
    palette[pmax(1, pmin(100, round(intensity[[i]] * 99) + 1))]
  }, character(1))
}

reorder_within <- function(x, by, within, fun = mean, sep = "___")
  stats::reorder(paste(x, within, sep = sep), by, FUN = fun)

scale_y_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  scale_y_discrete(..., labels = function(x) gsub(reg, "", x))
}

# %% Trait count and rate plot
#| Construct count bars, trait labels and lenient union-rate points.
fig4_plot_traits <- function(combined, label_top = Inf, label_char_frac = 0.0165,
                             label_gap = 0.03, label_size = 6 / (72.27 / 25.4)) {
  if (length(label_top) != 1 || is.na(label_top) || label_top < 0 ||
      (is.finite(label_top) && label_top != as.integer(label_top)))
    stop("label_top must be a non-negative integer, or Inf for every trait.",
         call. = FALSE)
  colors <- complete_category_colors(combined$trait_category)
  plotdat <- combined %>%
    group_by(trait_category) %>%
    mutate(
      intensity = if (max(count) == min(count)) rep(0.75, n())
                  else (count - min(count)) / (max(count) - min(count)),
      intensity = 0.25 + 0.75 * intensity,
      fill_col  = make_fill(trait_category, intensity, colors)) %>%
    ungroup() %>%
    mutate(
      trait_category = forcats::fct_drop(forcats::fct_reorder(trait_category, count, .fun = sum)),
      trait_raw      = trait,
      trait          = reorder_within(trait, count, trait_category),
      base_col       = unname(colors[as.character(trait_category)]))

  long <- bind_rows(
    plotdat %>% transmute(trait_category, trait, base_col, fill_col, metric = "a_count", value = count),
    plotdat %>% transmute(trait_category, trait, base_col, fill_col, metric = "b_rate",  value = pct_coloc_union))
  xline <- -max(plotdat$count) * 0.03
  seg  <- plotdat %>% count(trait_category, name = "nbar") %>% mutate(metric = "a_count")
  vref <- tibble(metric = "b_rate", xi = median(plotdat$pct_coloc_union, na.rm = TRUE))
  metric_labels <- c(a_count = "Colocalizing Genes", b_rate = "CS coloc rate ")

  count_max <- max(plotdat$count)
  count_scale <- scale_x_continuous()
  label_layer <- list()
  if (label_top > 0) {
    labdat <- plotdat %>%
      group_by(trait_category) %>%
      # slice_max() needs a constant n, so Inf ("label everything") becomes the
      # row count, which is an upper bound on any single category's size.
      slice_max(count, n = min(label_top, nrow(plotdat)), with_ties = FALSE) %>%
      ungroup() %>%
      mutate(metric = "a_count", label = trait_label(trait_raw))
    # The current print layout places all labels at count-axis position 2000.
    label_x <- 2000
    labdat$x <- label_x
    # Retain the four count-axis ticks used in the current figure.
    count_scale <- scale_x_continuous(breaks = c(0, 400, 800, 1200))
    # The harmonized coloc table carries near-duplicate trait strings from
    # different GWAS sources ("Asthma" and "asthma"), which collapse to the same
    # label on adjacent rows. That is what the data says, but it reads as a
    # mistake, so name them rather than letting them pass silently.
    dup_labels <- labdat %>% count(trait_category, label) %>% filter(n > 1)
    if (nrow(dup_labels) > 0)
      message("Repeated trait label(s) within a category: ",
              paste(sprintf("%s/%s x%d", dup_labels$trait_category,
                            dup_labels$label, dup_labels$n), collapse = ", "))
    label_layer <- geom_text(data = labdat, inherit.aes = FALSE,
                             aes(x = x, y = trait, label = label),
                             colour = "grey25", hjust = 1, size = label_size,
                             family = "Helvetica", show.legend = FALSE)
  }
  category_labeller <- c(GROUP_LABELS, setNames(
    setdiff(as.character(unique(plotdat$trait_category)), names(GROUP_LABELS)),
    setdiff(as.character(unique(plotdat$trait_category)), names(GROUP_LABELS))))

  ggplot(long, aes(y = trait)) +
    geom_col(data = ~ subset(.x, metric == "a_count"),
             aes(x = value, fill = fill_col), width = 0.85, show.legend = FALSE) +
    label_layer +
    geom_segment(data = seg, inherit.aes = FALSE,
                 aes(y = 0.5, yend = nbar + 0.5, x = xline, xend = xline),
                 color = "grey60", linewidth = 0.2) +
    geom_vline(data = vref, aes(xintercept = xi),
               linetype = "dashed", color = "grey70", linewidth = 0.15) +
    geom_segment(data = ~ subset(.x, metric == "b_rate"),
                 aes(x = 0, xend = value, yend = trait, color = base_col), linewidth = 0.15) +
    geom_point(data = ~ subset(.x, metric == "b_rate"),
               aes(x = value, color = base_col), size = 0.65) +
    scale_fill_identity() + scale_color_identity() + scale_y_reordered() +
    facet_grid(trait_category ~ metric, scales = "free", space = "free_y", switch = "y",
               labeller = labeller(
                 trait_category = as_labeller(category_labeller, default = label_value),
                 metric = as_labeller(metric_labels))) +
    ggh4x::facetted_pos_scales(x = list(
      count_scale,
      scale_x_continuous(labels = function(x) paste0(x, "%"),
                         expand = expansion(mult = c(0, 0.08))))) +
    ggh4x::force_panelsizes(cols = c(3, 1.6)) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_classic() + fig4_base_theme() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          panel.spacing.x = grid::unit(3, "pt"), panel.spacing.y = grid::unit(1.5, "pt"),
          strip.placement = "outside", strip.background = element_blank(),
          # Bold on the column headers only. The trait-category strips down the
          # side are already set apart by size and the grey rule beside them, and
          # bolding them too made the panel read as all heading.
          strip.text.y.left = element_text(angle = 0, hjust = 1, size = 6.5),
          strip.text.x = element_text(size = 7, face = "bold"),
          plot.margin = margin(2, 5, 2, 2, unit = "pt"))
}


# %% Rare-variant prediction plot
#| Use the original 30-bin histogram and the original percent-axis limits.
fig4_plot_r2 <- function(data) {
  ggplot(data, aes(x = 100 * difference)) +
    geom_histogram(colour = "white", linewidth = 0.14, bins = 30) +
    scale_x_continuous(limits = c(-20, 40)) +
    labs(x = expression(R^2 ~ increase ~ "(%)"), y = "Number of genes") +
    theme_classic() + fig4_base_theme() +
    theme(axis.title.x = element_text(size = 8.7))
}

# %% Conditional TWAS QQ plot
#| Compare calculated QQ coordinates. Null retains the current gray color.
fig4_plot_qq <- function(data) {
  ggplot(data, aes(x = expected, y = observed_display, colour = group)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.2, colour = "#666666") +
    geom_point(size = 0.4, alpha = 0.8) +
    facet_wrap(~test, scales = "free_y", labeller = as_labeller(c(common = "Common", rare = "Rare"))) +
    scale_colour_manual(values = c(Null = "grey50", `Colocalized genes` = "#275D8E"), name = NULL) +
    guides(colour = guide_legend(override.aes = list(shape = 16, size = 2.5, alpha = 1))) +
    labs(x = expression(Expected ~ TWAS ~ -log[10](P)), y = expression(Observed ~ TWAS ~ -log[10](P))) +
    theme_classic(base_size = 12) + fig4_base_theme() +
    theme(strip.background = element_blank(), panel.spacing.x = grid::unit(5, "pt"),
          panel.spacing.y = grid::unit(2, "lines"), axis.title = element_text(size = 8.7),
          legend.position = "none", legend.key.width = grid::unit(1.2, "lines"))
}

# %% Gene-set enrichment plot
#| Plot the fitted odds ratios with their confidence limits.
fig4_plot_enrichment <- function(data) {
  ggplot(data, aes(OR, reorder(type, OR), colour = term)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60", linewidth = 0.2) +
    geom_pointrange(aes(xmin = conf.low, xmax = conf.high), position = position_dodge(width = 0.55),
                    linewidth = 0.3, size = 0.65, stroke = 0.22) +
    scale_colour_manual(values = c(Rare = "#DB3A2FFF", Common = "grey60"), name = NULL) +
    labs(x = "Odds ratio per unit increase in Z-score", y = "") +
    theme_classic(base_size = 12) + fig4_base_theme() +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          strip.background = element_blank(), strip.text = element_text(face = "bold"),
          legend.position = c(0.9, 0.2), panel.spacing = grid::unit(1.4, "lines"),
          axis.title.x = element_text(size = 7), legend.key.size = grid::unit(3, "mm"))
}

# %% Figure 4 assembly
#| Keep the current panel layout. An unavailable panel A leaves an empty panel.
fig4_assemble <- function(a, b, c, d) {
  left <- if (is.null(a)) plot_spacer() else a + labs(tag = "A")
  top_right <- ((b + labs(tag = "B")) | (c + labs(tag = "C"))) + plot_layout(widths = c(0.85, 1.5))
  right <- wrap_elements(full = top_right) / (d + labs(tag = "D")) + plot_layout(heights = c(1, 0.75))
  (left | right) + plot_layout(widths = c(1.1, 1))
}
