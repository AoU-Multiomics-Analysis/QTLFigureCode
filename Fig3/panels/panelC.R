# Panel C — effect of -log(MAF) on cross-assay effect-size correlation:
# the interaction coefficient and 95% CI for each assay pair.
# Inputs: enloc_variant, annotations, susie_eqtl, susie_sqtl, susie_pqtl
source("fig3/00_setup.R")

# The original fit these models TWICE: once on the raw scale without an NA
# filter, then again with log1p() and NA filtering, discarding the first result.
# Only the second is kept.
#
# NOTE: the eQTL-pQTL model omits abs() on both posteriors, unlike the other two
# pairs. That asymmetry is preserved verbatim from the original; see NOTES.md,
# which flags it as possibly unintended (log1p() of a value below -1 is NaN).
FIT_FORMULA <- list(
  `eQTL-sQTL` = log1p(abs(posterior_mean.y)) ~ log1p(abs(posterior_mean.x)) * log_maf,
  `eQTL-pQTL` = log1p(posterior_mean.y)      ~ log1p(posterior_mean.x)      * log_maf,
  `sQTL-pQTL` = log1p(abs(posterior_mean.y)) ~ log1p(abs(posterior_mean.x)) * log_maf
)

fit_pair <- function(d, pair) {
  d <- filter(d, !is.na(posterior_mean.x), !is.na(posterior_mean.y))
  broom::tidy(lm(FIT_FORMULA[[pair]], data = d), conf.int = TRUE) %>%
    filter(!str_detect(term, "Int")) %>%
    mutate(type = pair)
}

interaction_res <- purrr::imap_dfr(coloc_model_input(), ~ fit_pair(.x, .y))

pC <- interaction_res %>%
  filter(str_detect(term, ":")) %>%
  ggplot(aes(x = estimate, y = type)) +
    geom_linerange(aes(xmin = conf.low, xmax = conf.high), color = "grey40") +
    geom_point(size = 4, color = "#364B63") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    labs(x = "Effect of -log(MAF) on cross-ome correlation", y = NULL) +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          plot.margin = margin(l = 60, unit = "pt"))

saveRDS(pC, file.path(OBJ_DIR, "C.rds"))
message("Panel C built.")
