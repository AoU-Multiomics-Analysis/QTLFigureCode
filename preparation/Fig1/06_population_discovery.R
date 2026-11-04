# Recompute S4 proportions from the 21 source permutation-result tables.
# This historical panel uses BH FDR per file. Main Figure 1 novelty uses source qval.
source("preparation/common.R")
config <- prep_init()
assays <- c("eQTL", "sQTL", "pQTL")
groups <- c("COMB", "AFR", "AMR", "EAS", "EUR", "SAS", "MID")
target <- "Supplemental Figures/01_Supplemental_drafts/Figure1_QTL_discovery_by_population/population_QTL_discovery_data.tsv"
samples <- read_tsv(prep_input("input", "reference/COMB_eQTL_9k_sample_list.tsv"),
  col_types = cols_only(research_id = col_character()), progress = FALSE)
ancestry <- read_tsv(prep_input("input", "echo_v4_r2.ancestry_preds.tsv"),
  col_types = cols_only(research_id = col_character(), ancestry_pred_other = col_character()), progress = FALSE) %>%
  semi_join(samples, by = "research_id")
stopifnot(nrow(samples) == 8895L, !anyDuplicated(samples$research_id),
          nrow(ancestry) == nrow(samples), !anyDuplicated(ancestry$research_id), !anyNA(ancestry))
reference_n <- ancestry %>% mutate(population = str_to_upper(ancestry_pred_other)) %>% count(population) %>%
  bind_rows(tibble(population = "COMB", n = nrow(samples)))
fills <- c(COMB = "#CF597E", AFR = "#C8D892", AMR = "#7FC385", EAS = "#ECA676",
           EUR = "#30AC86", SAS = "#EBCE8D", MID = "#E47D72")
summaries <- crossing(assay = assays, population = groups) %>% pmap_dfr(function(assay, population) {
  relative <- if (assay == "pQTL") paste0("Fig1Data/pQTL/AoU_10750_", population, "_pQTL.cis_qtl.txt.gz") else
    paste0("Fig1Data/", assay, "/", population, ".cis_qtl.txt.gz")
  message("S4: Calculate BH FDR within ", assay, " / ", population)
  d <- read_tsv(prep_input("input", relative),
    col_types = cols_only(pval_beta = col_character()), progress = FALSE)
  fdr <- p.adjust(as.double(d$pval_beta), method = "BH")
  stopifnot(!anyNA(fdr), all(between(fdr, 0, 1)))
  tibble(assay, population, proportion_significant = mean(fdr < .05), source_n_tested = nrow(d),
         significant_count_from_saved_summary = sum(fdr < .05))
}) %>% left_join(reference_n, by = "population", relationship = "many-to-one") %>%
  transmute(assay = factor(assay, levels = assays),
    population = factor(recode(population, COMB = "Combined"), levels = c("Combined", groups[-1])),
    proportion_significant, source_n_tested, significant_count_from_saved_summary,
    source_reference_sample_n = n, colour = unname(fills[as.character(recode(population, Combined = "COMB"))]),
    marker = if_else(population == "Combined", 23L, 21L)) %>% arrange(assay, population)
saved <- read_tsv(prep_reference("paper", target), show_col_types = FALSE) %>%
  mutate(assay = factor(assay, levels = assays),
    population = factor(population, levels = c("Combined", groups[-1]))) %>% arrange(assay, population)
stopifnot(isTRUE(all.equal(summaries, saved, tolerance = 1e-12, check.attributes = FALSE)))
prep_write_tsv(summaries, "paper", target)
prep_write_tsv(summaries, "paper", "output/figure1_population_qtl_discovery/population_QTL_discovery_data.tsv")
message("S4 verified: all 21 raw-result proportions, tested-row counts, significant counts and reference sample counts match.")
