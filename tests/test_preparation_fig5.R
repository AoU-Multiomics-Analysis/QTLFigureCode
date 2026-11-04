suppressPackageStartupMessages(library(tidyverse))
source("preparation/Fig5/functions.R")

reference <- tibble(
  targetId = c("ENSG1", "ENSG1", "ENSG2", "ENSG3"),
  targetFromSourceId = c("A", "A", "B", "C"),
  confidence = c("Strong", "Strong", "Definitive", "Limited"),
  allelicRequirements = list("AD", "AR", "AR", "AR")
)
ref <- fig5_reference(reference)
raw <- tibble(
  variant = c("v1", "v1", "v2", "v3", "v4", "v5"),
  molecular_trait_id = c("ENSG1.1", "ENSG1.1", "ENSG2.1", "ENSG2.1",
    "ENSG3.1", "ENSG2.1"),
  assay = c("eQTL", "pQTL", "eQTL", "sQTL", "eQTL", "eQTL"),
  gvs_max_af = c(.009, .009, .01, 0, .001, .001),
  pip = c(.95, .98, 1, .9, 1, .99),
  posterior_mean = c(-.4, .5, -.5, -.3, -.8, .5),
  clinvar_classification = c("", "", "benign", "uncertain significance", "", "pathogenic")
)
pairs <- fig5_select_pairs(raw, ref)
stopifnot(setequal(pairs$variant, c("v1", "v3", "v5")))
frequencies <- tibble(
  variant = c("v1", "v3", "v5"),
  gvs_afr_af = c(.008, .004, .004), gvs_amr_af = c(.008, .004, .001),
  gvs_eas_af = c(.001, .001, .001), gvs_eur_af = c(.006, .002999, .003),
  gvs_mid_af = 0, gvs_sas_af = 0
)
af <- fig5_assign_populations(frequencies)
stopifnot(
  af$population_A[af$variant == "v1"] == "afr",
  af$tied_populations[af$variant == "v1"] == "afr;amr",
  af$population_B[af$variant == "v1"] == "eur",
  af$population_B[af$variant == "v3"] == "afr",
  af$population_B[af$variant == "v5"] == "eur"
)
assignments <- fig5_variant_assignments(pairs, frequencies, ref)
summaries <- fig5_counts_and_classes(assignments, pairs)
stopifnot(
  nrow(assignments) == 4L,
  sum(summaries$counts$total[summaries$counts$allelicRequirements == "AD"]) == 1,
  sum(summaries$counts$total[summaries$counts$allelicRequirements == "AR"]) == 3,
  sum(summaries$proportions$n) == 4,
  fig5_clinvar_class(c(NA, "", "benign,uncertain significance,likely pathogenic",
    "benign,uncertain significance", "likely benign")) %>%
    identical(c("Unannotated", "Unannotated", "P/LP", "VUS", "B/LB"))
)
bad_af <- frequencies
bad_af$gvs_afr_af[[1]] <- NA_real_
stopifnot(inherits(try(fig5_assign_populations(bad_af), silent = TRUE), "try-error"))

funnel_raw <- tibble(
  variant = c("reference", "candidate", "boundary", "lowpip", "reference"),
  molecular_trait_id = "ENSG2.1", assay = c("eQTL", "eQTL", "eQTL", "eQTL", "pQTL"),
  gvs_max_af = .001, pip = c(.95, .96, .99, .9, .99),
  posterior_mean = c(-.4, .5, -.3, -.9, .2),
  clinvar_classification = c("pathogenic", "uncertain significance", "", "", "pathogenic")
)
funnel <- fig5_funnel(funnel_raw, ref)
stopifnot(
  identical(funnel$stages$variants, c(4L, 3L, 1L, 3L)),
  identical(funnel$stages$pathogenic_variants, c(1L, 1L, 1L, 1L)),
  funnel$stages$genes[[4]] == 1L,
  "candidate" %in% funnel$selected_pairs$variant,
  any(funnel$selected_pairs$posterior_mean > 0),
  is.na(funnel$stages$enrichment[[4]])
)
membership <- tibble(variant = c("r", "r", "a", "b", "c"),
  is_pathogenic = c(TRUE, TRUE, FALSE, FALSE, FALSE))
test <- fig5_subset_test(membership, tibble(variant = c("r", "r", "a")),
  "selected", "universe", "fixture")
stopifnot(test$retained_pathogenic == 1L, test$retained_other == 1L,
  test$not_retained_pathogenic == 0L, test$not_retained_other == 2L,
  test$universe_variants == 4L)
thresholds <- fig5_threshold_tables(tibble(
  variant = c("one", "one", "two", "three"),
  gene = c("a", "b", "c", "d"), log2_aFC = c(-1, -2, 0, -.5)
), c(0, 50, 75, 100))
stopifnot(identical(thresholds$counts$variants, c(2L, 1L, 0L, 0L)),
  all(thresholds$counts$total_assessed == 3L))
locus_source <- tibble(variant = c("a", "b"), position = c(51912316, 51914316),
  pip = c(.9, .95), cs_id = c("set1", "set2"), gvs_max_af = c(.01, .009))
locus_template <- locus_source %>%
  mutate(cs_id = factor(cs_id), rare = c(FALSE, TRUE), causal = TRUE,
    in_zoom = c(TRUE, FALSE), type = "Expression", annotation = c("kept", "kept"))
locus_plot <- ggplot(locus_template, aes(position, pip)) +
  geom_point(data = locus_template %>% slice(1))
refreshed <- fig5_refresh_locus_template(locus_plot, locus_source %>% slice(2, 1))
stopifnot(identical(refreshed$data$variant, locus_template$variant),
  identical(refreshed$data$annotation, locus_template$annotation),
  all(refreshed$data$causal), !refreshed$data$rare[[1]],
  refreshed$data$rare[[2]], nrow(refreshed$layers[[1]]$data) == 1L)
missing <- locus_source %>% slice(1)
changed <- locus_source %>% mutate(pip = pip - .01)
extra <- bind_rows(locus_source, locus_source %>% slice(1) %>% mutate(variant = "new"))
stopifnot(inherits(try(fig5_refresh_locus_template(locus_plot, missing), silent = TRUE), "try-error"),
  inherits(try(fig5_refresh_locus_template(locus_plot, changed), silent = TRUE), "try-error"),
  inherits(try(fig5_refresh_locus_template(locus_plot, extra), silent = TRUE), "try-error"))
cat("Figure 5 preparation tests passed.\n")
