# =============================================================================
# Fine-mapping annotation enrichment, modelled as logistic regression on PIP.
#
# This is an ALTERNATIVE to the enrichment behind Figure 1 panel H (rendered
# tag I), not a replacement, and it writes its own standalone figure. The two
# ask different questions:
#
#   panel H  Fisher exact test on a 2x2 table of (PIP > 0.9) x (annotated),
#            after throwing away every variant with 0.01 <= PIP <= 0.9. The
#            effect is an odds ratio: how much more often a putatively causal
#            variant carries the annotation than a background variant does.
#
#   here     glm(annotated ~ PIP, binomial) over EVERY variant in the table.
#            The effect is a log-odds slope: how much the log odds of carrying
#            the annotation rises as PIP goes from 0 to 1. No variant is
#            discarded and the middle of the PIP range informs the fit.
#
# A positive slope and an odds ratio above one mean the same thing
# qualitatively, but the numbers are on different scales and are NOT
# interchangeable.
#
# Run from the project root, with the SuSiE tables staged:
#   FIG1_INPUT=/path/to/susie_tables Rscript fig1/enrichment_glm.R
# =============================================================================
source("fig1/00_setup.R")

# --- the model ---------------------------------------------------------------
# Shared with panel H, and memoized, so whichever runs first pays for the fit.
source("fig1/enrichment_glm_fit.R")

results <- FitAnnotationEnrichmentGlm()

sep <- results %>% filter(separated)
if (nrow(sep)) {
  message("Excluded from the plot -- separated fits with no finite slope:\n",
          paste0("  ", sep$type, " ", sep$name,
                 "  estimate=", signif(sep$estimate, 3),
                 "  se=", signif(sep$std.error, 3), collapse = "\n"))
}

# --- the figure --------------------------------------------------------------
source("fig1/enrichment_glm_plot.R")
