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

# Fits every assay and returns one row per annotation x assay slope. Memoized to
# fig1/objects/enrichment_glm.tsv, because both panel H and the standalone
# fig1/enrichment_glm.R figure want the same fits and the eQTL table takes a few
# minutes. Pass refresh = TRUE after changing the model or the inputs.
FitAnnotationEnrichmentGlm <- function(refresh = FALSE) {
  cache <- file.path(OBJ_DIR, "enrichment_glm.tsv")
  if (!refresh && file.exists(cache)) {
    message("Reusing cached fits: ", cache, "  (refresh = TRUE to refit)")
    return(data.table::fread(cache))
  }

  annot <- input_files("annot")
  fits <- purrr::pmap_dfr(
    list(annot$path, annot$Assay),
    function(p, ty) {
      message("Fitting ", ty, " (", basename(p), ")")
      GlmAnnotationEnrichment(p) %>% mutate(type = ty)
    }
  )

  # FDR across every annotation x assay slope, matching how panel H adjusts.
  results <- fits %>%
    filter(term == "pip") %>%
    ungroup() %>%
    mutate(
      padj      = p.adjust(p.value, method = "fdr"),
      conf.low  = estimate - 1.96 * std.error,
      conf.high = estimate + 1.96 * std.error,
      # A perfectly separated fit (every annotated variant above some PIP, or
      # none at all) has no finite maximum-likelihood slope. glm still returns
      # one, but it is an artefact of where the optimiser stopped and its
      # standard error is enormous.
      separated = std.error > 10
    ) %>%
    arrange(desc(estimate))

  data.table::fwrite(results, cache, sep = "\t")
  message("Wrote ", cache, " (", nrow(results), " annotation x assay slopes)")
  results
}
