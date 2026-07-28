# =============================================================================
# VEP consequence expansion, shared by Fig1 panel H and Fig3 panel E.
#
# VEP reports every consequence for a variant as ONE comma-separated string
# ("intron_variant,NMD_transcript_variant"), and ~20% of rows carry more than
# one term. Enrichment code picks annotation columns up via where(is.logical),
# so unless each term becomes its own logical column the whole coding / splice
# / UTR / noncoding-RNA half of the annotation space is silently absent.
#
# One-hot expansion is done by matrix construction rather than separate_rows()
# + pivot_wider(): these tables carry ~80 columns, and pivot_wider would treat
# every one of them as an id column.
# =============================================================================
expand_vep_consequences <- function(d) {
  if (!"consequence" %in% names(d) || !nrow(d)) return(d)
  terms <- lapply(strsplit(as.character(d$consequence), ",", fixed = TRUE),
                  function(x) unique(trimws(x)))
  lvl <- setdiff(unique(unlist(terms)), c("", ".", NA))
  if (!length(lvl)) return(d)
  # Never clobber a column that already exists (the cCRE flags are already
  # one-hot); a duplicate name would make select() ambiguous downstream.
  lvl <- setdiff(lvl, names(d))
  if (!length(lvl)) return(d)
  m <- vapply(lvl, function(t) vapply(terms, function(x) t %in% x, logical(1)),
              logical(length(terms)))
  dplyr::bind_cols(d, as.data.frame(m, stringsAsFactors = FALSE))
}
