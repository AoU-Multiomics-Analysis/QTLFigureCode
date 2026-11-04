# %% Count reference samples and test discovery in each source table
#| The historical panel uses BH FDR within each permutation-result file. Main Figure 1 novelty uses source qval.
#| The sample sizes are common RNA reference-cohort counts reused across assays.
assays <- c("eQTL", "sQTL", "pQTL")
groups <- c("COMB", "AFR", "AMR", "EAS", "EUR", "SAS", "MID")
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
s4 <- crossing(assay = assays, population = groups) %>% pmap_dfr(function(assay, population) {
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

message("Calculated discovery proportions and counts from 21 source result tables.")

# %% Draw the calculated discovery proportions
#| Each point is a population and assay summary. The x axis uses the counted reference sample count. It is not an assay-specific participant count. No confidence interval or fitted curve is added.
suppressPackageStartupMessages(library(ggrepel))
fills <- s4 %>% distinct(population, colour) %>% deframe()
shapes <- s4 %>% distinct(population, marker) %>% deframe()
figure <- ggplot(s4, aes(source_reference_sample_n, proportion_significant)) +
  geom_point(aes(fill = population, shape = population), size = 2.3, colour = "#263B45", stroke = 0.35) +
  geom_text_repel(aes(label = population), seed = 71, size = 2.6,
                  min.segment.length = 0, box.padding = 0.3,
                  point.padding = 0.2, max.overlaps = Inf) +
  facet_grid(~assay) +
  scale_fill_manual(values = fills) +
  scale_shape_manual(values = shapes) +
  scale_x_log10(limits = c(280, 11000), breaks = c(300, 1000, 3000, 10000), labels = scales::label_comma()) +
  scale_y_continuous(limits = c(0, 1), labels = scales::label_percent(), breaks = seq(0, 1, 0.25)) +
  labs(x = "Reference cohort sample size (log scale)", y = "Tested genes with significant QTLs (%)") +
  theme_classic(base_size = 9) +
  theme(strip.background = element_blank(), legend.position = "none",
        panel.spacing.x = grid::unit(5, "mm"))
print(figure)

# %% Save aggregate outputs
#| Save the figure and the 21 study summaries. No participant records are exported. The source field source_n_tested counts tested rows, not significant genes.
save_plot(figure, "supplements/S4", 184, 85)
readr::write_tsv(s4 %>% select(assay, population, proportion_significant, source_n_tested,
                              significant_count_from_saved_summary, source_reference_sample_n),
                 output_file("supplements/S4_summary.tsv"))
