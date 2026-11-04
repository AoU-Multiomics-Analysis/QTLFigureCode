suppressPackageStartupMessages(library(tidyverse))
stopifnot(file.exists("notebooks/current/analysis/S15_helpers.R"))
source("notebooks/current/analysis/S15_helpers.R")
notebook_paths <- c("notebooks/current/supplements/S15.R", "notebooks/current/analysis/S15_helpers.R")
notebook_code <- paste(unlist(lapply(notebook_paths, readLines, warn = FALSE)), collapse = "\n")
stopifnot(!grepl("readRDS|saveRDS|prep_output|prep_save_rds|parse\\(|eval\\(|source\\(", notebook_code))

# A missing VEP consequence must remain in the background for every flag fit.
fixture <- tibble(variant = paste0("v", 1:16), pip = rep(c(0, .2, .7, 1), each = 4),
  consequence = rep(c("missense_variant", NA, ".", "3_prime_UTR_variant"), 4),
  CTCF.only = rep(c(TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE), 2),
  CTCF.bound = rep(c(TRUE, FALSE), 8), pELS = rep(c(FALSE, TRUE), 8),
  DNase.H3K4me3 = rep(c(TRUE, TRUE, FALSE, FALSE), 4),
  dELS = rep(c(FALSE, FALSE, TRUE, TRUE), 4), PLS = rep(c(TRUE, FALSE, TRUE, FALSE), 4),
  FANTOM5 = rep(c(FALSE, TRUE, FALSE, TRUE), 4))
path <- tempfile(fileext = ".tsv")
readr::write_tsv(fixture, path)
fits <- GlmAnnotationEnrichment(path)
expected <- broom::tidy(glm(CTCF.only ~ pip, family = binomial(), data = fixture))
actual <- fits %>% filter(name == "CTCF.only") %>% arrange(term)
stopifnot(isTRUE(all.equal(actual$estimate, expected %>% arrange(term) %>% pull(estimate), tolerance = 1e-12)),
  !NO_CONSEQUENCE %in% fits$name, !anyNA(fits$estimate))
summary <- s15_model_summaries(bind_rows(mutate(fits, type = "eQTL"), mutate(fits, type = "sQTL")))
stopifnot(all(summary$term == "pip"),
  isTRUE(all.equal(summary$conf.low, summary$estimate - 1.96 * summary$std.error)),
  isTRUE(all.equal(summary$padj, p.adjust(summary$p.value, method = "fdr"))))
plot_data <- s15_plot_data(summary)
stopifnot(all(plot_data$name %in% names(PANEL_ANNOTATION_LABELS)),
  !anyNA(plot_data$functional_group), all(plot_data$lower > 0),
  all(plot_data$upper >= plot_data$lower))
unlink(path)
message("Inline S15 annotation-background and model-summary checks passed.")
