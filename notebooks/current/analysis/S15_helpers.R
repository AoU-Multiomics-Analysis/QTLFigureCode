# =============================================================================
# The model behind fig1/enrichment_glm.R. Defines only -- sourcing this runs
# nothing -- so tests can exercise the fit without the controlled-access tables.
# =============================================================================

# The genuine logical annotation flags, named rather than picked up with
# `where(is.logical)`. They cannot be found by type: in the SPLICING table tss,
# gene_name, distTSS, Constrained and lof.pLI are NA for every single row, so
# fread types all five as logical and `where(is.logical)` sweeps them in as
# annotations. One of them is then all-NA after the pivot, which leaves glm()
# no response to fit and aborts the run. Naming the columns also means a renamed
# or missing flag is an error here instead of an annotation that quietly
# disappears from the results, which is how two UTR annotations once went
# missing from panel H (see NOTES.md).
CCRE_FLAG_COLS <- c("CTCF.only", "CTCF.bound", "pELS", "DNase.H3K4me3",
                    "dELS", "PLS", "FANTOM5")

# Where variants with no VEP consequence are parked. They must stay in the
# table -- a variant with no consequence is still a background variant for
# every cCRE flag, and dropping the rows would quietly shrink the denominator
# of every other fit -- but "no consequence recorded" is not an annotation to
# regress, so the column is removed once the one-hot pivot has happened.
NO_CONSEQUENCE <- ".no_consequence_recorded"

# One row per (variant, PIP, annotation): the cCRE flags arrive as logical
# columns and VEP consequences as one comma-separated string, so the
# consequences are one-hot expanded and then folded in with the flags. Returns
# broom::tidy() output per annotation, i.e. an (Intercept) row and a pip row.
GlmAnnotationEnrichment <- function(path) {
  d <- data.table::fread(path,
                         select = c("variant", "pip", "consequence", CCRE_FLAG_COLS))

  missing <- setdiff(CCRE_FLAG_COLS, names(d))
  if (length(missing))
    stop("annotation flags absent from ", basename(path), ": ",
         paste(missing, collapse = ", "))
  not_logical <- CCRE_FLAG_COLS[!vapply(d[, ..CCRE_FLAG_COLS], is.logical, logical(1))]
  if (length(not_logical))
    stop("annotation flags are not logical in ", basename(path), ": ",
         paste(not_logical, collapse = ", "))

  message("  ", nrow(d), " rows, ", dplyr::n_distinct(d$variant), " distinct variants")

  d %>%
    select(variant, pip, where(is.logical), consequence) %>%
    separate_rows(consequence, sep = ",") %>%
    # NA, "" and VEP's "." placeholder all mean the same thing. The empty string
    # matters in its own right: pivot_wider cannot create a column named "" and
    # errors out, whereas NA silently becomes a column called "NA".
    mutate(consequence = ifelse(is.na(consequence) | consequence %in% c("", "."),
                                NO_CONSEQUENCE, consequence)) %>%
    mutate(value = TRUE) %>%
    distinct() %>%
    pivot_wider(names_from = consequence, values_from = value, values_fill = FALSE) %>%
    select(-any_of(NO_CONSEQUENCE)) %>%
    pivot_longer(where(is.logical)) %>%
    group_by(name) %>%
    # An annotation that is TRUE for every variant, or FALSE for every variant,
    # has nothing to regress; glm() errors rather than returning a null fit.
    filter(dplyr::n_distinct(value) == 2) %>%
    do(broom::tidy(glm(value ~ pip, family = "binomial", data = .))) %>%
    arrange(desc(estimate))
}

# Adjust all annotation-by-assay PIP slopes together, before the plot subset.
s15_model_summaries <- function(fits) {
  fits %>% filter(term == "pip") %>% ungroup() %>%
    mutate(padj = p.adjust(p.value, method = "fdr"),
      conf.low = estimate - 1.96 * std.error,
      conf.high = estimate + 1.96 * std.error,
      separated = std.error > 10) %>% arrange(desc(estimate))
}

# These labels select the same twelve annotations as Figure 1H.
PANEL_ANNOTATION_LABELS <- c(
  splice_donor_variant = "Splice donor", splice_acceptor_variant = "Splice acceptor",
  stop_gained = "Stop gained", missense_variant = "Missense",
  frameshift_variant = "Frameshift", start_lost = "Start lost",
  PLS = "Promoter", dELS = "Distal enhancer", `5_prime_UTR_variant` = "5' UTR",
  `3_prime_UTR_variant` = "3' UTR", pELS = "Proximal enhancer",
  splice_region_variant = "Splice region")

variant_functional_groups <- tibble::tribble(
  ~name, ~functional_group,
  "3_prime_UTR_variant", "Regulatory Annotation",
  "5_prime_UTR_variant", "Regulatory Annotation",
  "CTCF.bound", "Regulatory Annotation",
  "DNase.H3K4me3", "Regulatory Annotation",
  "FANTOM5", "Regulatory Annotation",
  "PLS", "Regulatory Annotation", "pELS", "Regulatory Annotation", "dELS", "Regulatory Annotation",
  "coding_sequence_variant", "Coding Annotation", "frameshift_variant", "Coding Annotation",
  "inframe_deletion", "Coding Annotation", "inframe_insertion", "Coding Annotation",
  "missense_variant", "Coding Annotation", "synonymous_variant", "Coding Annotation",
  "start_lost", "Coding Annotation", "stop_gained", "Coding Annotation",
  "stop_lost", "Coding Annotation", "stop_retained_variant", "Coding Annotation",
  "splice_acceptor_variant", "Splice Annotation", "splice_donor_variant", "Splice Annotation",
  "splice_donor_5th_base_variant", "Splice Annotation", "splice_donor_region_variant", "Splice Annotation",
  "splice_region_variant", "Splice Annotation", "splice_polypyrimidine_tract_variant", "Splice Annotation",
  "non_coding_transcript_exon_variant", "Noncoding RNA", "non_coding_transcript_variant", "Noncoding RNA",
  "mature_miRNA_variant", "Noncoding RNA", "downstream_gene_variant", "Intronic/Intergenic",
  "upstream_gene_variant", "Intronic/Intergenic", "intron_variant", "Intronic/Intergenic",
  "intergenic_variant", "Intronic/Intergenic")

s15_plot_data <- function(results) {
  data <- results %>% filter(name %in% names(PANEL_ANNOTATION_LABELS)) %>%
    left_join(variant_functional_groups, by = "name", relationship = "many-to-one") %>%
    mutate(annotation_clean = unname(PANEL_ANNOTATION_LABELS[name])) %>%
    group_by(functional_group) %>% mutate(row = as.numeric(factor(annotation_clean))) %>% ungroup()
  stopifnot(!anyNA(data$functional_group))
  if (any(data$separated)) {
    message("Exclude selected annotation fits with no finite maximum-likelihood slope.")
    data <- filter(data, !separated)
  }
  data %>% mutate(odds_ratio = exp(estimate), lower = exp(conf.low), upper = exp(conf.high),
    type = factor(type, c("eQTL", "sQTL", "pQTL")), significant = padj < 0.05)
}
