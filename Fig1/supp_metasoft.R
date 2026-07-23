# =============================================================================
# Supplemental figure — cross-ancestry meta-analysis (MetaSoft) summary.
# NOT part of Figure 1; produces xGeneMetasoftSummary.pdf.
#
# Reshapes per-ancestry QTL permutation betas into MetaSoft input, runs MetaSoft
# (m-value mode) per assay, and summarises the resulting M-values two ways:
#   - proportion of genes by number of ancestry groups supporting the effect
#   - M-value distribution per population, faceted by assay
#
# EXTERNAL DEPENDENCY: MetaSoft (Java). Set METASOFT_JAR to the Metasoft.jar path.
#   Download: http://genetics.cs.ucla.edu/meta/repository/2.0.1/Metasoft.zip
# If the per-assay <ASSAY>.metasoft.tsv outputs already exist in OUTDIR, the
# MetaSoft run is skipped and the existing outputs are used.
#
# Run from project root:  Rscript fig1/supp_metasoft.R
# =============================================================================
source("fig1/00_setup.R")

METASOFT_JAR <- Sys.getenv("METASOFT_JAR", "Metasoft.jar")

# --- MetaSoft functions (folded from QTLAnalysisUtils, bug-fixed) -----------
FormatMetaSoft <- function(QTL_input) {
  QTL_input %>%
    filter(group != "COMB") %>%
    select(phenotype_id, slope, slope_se, group) %>%
    dplyr::rename("beta" = 2, "se" = 3) %>%
    pivot_wider(names_from = group, values_from = c(beta, se)) %>%
    na.omit() %>%
    dplyr::rename("SNP" = 1) %>%
    select(SNP, contains("AFR"), contains("SAS"), contains("EUR"),
           contains("AMR"), contains("EAS"), contains("MID"))
}

PlotMetaSoftMvaluesProps <- function(merged_metasoft) {
  merged_metasoft %>%
    select(RSID, V23, V24, V25, V26, V27, V28, type) %>%
    pivot_longer(!c(RSID, type)) %>%
    filter(value >= .9) %>%
    group_by(RSID, type) %>%
    summarize(mval_count = dplyr::n()) %>%
    group_by(type) %>%
    dplyr::count(mval_count) %>%
    mutate(prop = n / sum(n)) %>%
    ggplot(aes(x = mval_count, y = prop, fill = type)) +
    geom_col(position = "dodge") +
    xlab("Ancestry groups supporting effect  \n(Number of M vals >= .9 per gene)") +
    ylab("Proportion of genes") +
    paletteer::scale_fill_paletteer_d("ltc::minou")
}

# FIX: operate on the argument instead of a global `metasoft_all`.
PlotMetaSoftMvalueByPop <- function(merged_metasoft) {
  merged_metasoft %>%
    select(RSID, V23, V24, V25, V26, V27, V28, type) %>%
    na.omit() %>%
    dplyr::rename("AFR" = 2, "SAS" = 3, "EUR" = 4, "AMR" = 5, "EAS" = 6, "MID" = 7) %>%
    pivot_longer(!c(RSID, type)) %>%
    ggplot(aes(x = name, y = value, fill = name)) +
    geom_boxplot(show.legend = FALSE) +
    scale_fill_carto_d(palette = "Temps") +
    xlab("Population") +
    ylab("M-values") +
    facet_wrap(~type, nrow = 1)
}

# --- Read per-ancestry permutations via the manifest ------------------------
# The "perm" key covers all three assays; Assay is attached from the manifest.
assay_labels <- c(eQTL = "Expression", pQTL = "Protein", sQTL = "Splicing")
QTL_perm <- load_input("perm")

message("Preparing MetaSoft inputs ...")
for (a in names(assay_labels)) {
  QTL_perm %>%
    filter(Assay == assay_labels[[a]]) %>%
    FormatMetaSoft() %>%
    write_tsv(file.path(OUTDIR, paste0(a, "_wide_metasoft_input.tsv")), col_names = FALSE)
}

# --- Run MetaSoft (skipped if outputs already present) ----------------------
run_metasoft <- function(a) {
  out <- file.path(OUTDIR, paste0(a, ".metasoft.tsv"))
  if (file.exists(out)) { message("  using existing ", out); return(invisible()) }
  if (!file.exists(METASOFT_JAR))
    stop("MetaSoft jar not found at '", METASOFT_JAR, "'. Set METASOFT_JAR or pre-stage ",
         a, ".metasoft.tsv in OUTDIR.")
  system(paste("java -jar", shQuote(METASOFT_JAR), "-mvalue -input",
               file.path(OUTDIR, paste0(a, "_wide_metasoft_input.tsv")),
               "-output", out))
}
invisible(lapply(names(assay_labels), run_metasoft))

# --- Summarise and plot -----------------------------------------------------
metasoft_all <- purrr::map_dfr(names(assay_labels), function(a)
  fread(file.path(OUTDIR, paste0(a, ".metasoft.tsv"))) %>% mutate(type = a))

metasoft_plot_prop <- PlotMetaSoftMvaluesProps(metasoft_all)
metasoft_plot_pop  <- PlotMetaSoftMvalueByPop(metasoft_all)

# FIX: original composed the undefined `mvalue_per_pop` (the function's internal
# return-variable name); use the built plot object instead.
metasoft_summary <- (metasoft_plot_prop / metasoft_plot_pop) +
  plot_annotation(tag_levels = "A")

ggsave(file.path(OUTDIR, "xGeneMetasoftSummary.pdf"), metasoft_summary,
       width = 12, height = 7, dpi = DPI)
message("Wrote ", file.path(OUTDIR, "xGeneMetasoftSummary.pdf"))
