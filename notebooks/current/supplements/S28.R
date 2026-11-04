# %% Source selection and threshold functions
# @include preparation/Fig5/functions.R

# %% Select rare eQTL variants and calculate expression-reduction thresholds
#| Keep PIP >0.9, maximum AF <1%, ClinGen targets, and uncertain or empty ClinVar records. Deduplicate variants after choosing the strongest decreasing effect.
message("Read the original S28 expression and ClinGen inputs.")
raw <- as_tibble(data.table::fread(prep_input("input",
  "COMB.expression_SusieMerged.annotated.ancestry_skew.tsv.gz"),
  select = c("variant", "molecular_trait_id", "gvs_max_af", "pip",
    "clinvar_classification", "log2_aFC"), na.strings = "NA"))
reference <- arrow::read_parquet(prep_input("input",
  "part-00000-006c7897-e219-4036-abbf-af4482a5177e-c000.snappy.parquet")) %>%
  as_tibble() %>% fig5_reference()
selected <- raw %>%
  mutate(clinvar_classification = replace_na(clinvar_classification, ""),
    gene = str_extract(molecular_trait_id, "ENSG[0-9]+")) %>%
  filter(gvs_max_af < 0.01, pip > 0.9,
    str_detect(clinvar_classification, "uncertain") | clinvar_classification == "",
    gene %in% reference$targetId) %>%
  select(variant, gene, log2_aFC, clinvar_classification, pip, gvs_max_af)
stopifnot(nrow(selected) == 576L, n_distinct(selected$variant) == 574L,
  n_distinct(selected$gene) == 453L, all(is.finite(selected$log2_aFC)))
curve_result <- fig5_threshold_tables(selected, seq(0, 100, by = 1))
landmark_result <- fig5_threshold_tables(selected,
  sort(unique(c(seq(0, 100, by = 5), 100 * (1 - 2^-0.5), 99))))
curve <- curve_result$counts
landmarks <- landmark_result$counts
dots <- landmarks %>% filter(reduction_threshold_pct %in% c(0, 10, 20, 25, 40, 50, 60, 70, 75, 80, 90, 95, 100) |
                              abs(log2_afc_cutoff + 0.5) < 1e-10)
half <- landmarks %>% filter(reduction_threshold_pct == 50)
legacy <- landmarks %>% filter(abs(log2_afc_cutoff + 0.5) < 1e-10)
stopifnot(all(curve$total_assessed == 574L), all(diff(curve$variants) <= 0),
          half$variants == 45L, legacy$variants == 105L)

# %% Draw the calculated thresholds
#| Expression reduction is 100 × (1 − 2^log2_aFC). A variant passes a threshold when at least one qualifying target has a strictly greater reduction. Straight lines join selected threshold points as in the current source figure; they do not estimate the intervening empirical counts.
top_breaks <- 100 * (1 - 2^c(0, -0.5, -1, -2, -log2(10)))
figure <- ggplot(dots, aes(reduction_threshold_pct, variants)) +
  geom_vline(xintercept = c(legacy$reduction_threshold_pct, 50), colour = "#A8AFB4", linewidth = 0.3, linetype = "dashed") +
  geom_line(linewidth = 0.55, colour = "#145D69") +
  geom_point(size = 1.65, colour = "#145D69") +
  geom_point(data = bind_rows(legacy, half), size = 2.5, colour = "#D55E00") +
  annotate("text", x = 1, y = 270, label = "254", hjust = 0, size = 2.8, colour = "#145D69") +
  annotate("text", x = 98, y = 270, label = "574 variants assessed", hjust = 1, size = 2.8) +
  annotate("text", x = 32, y = 137, label = "105 variants (18.3%)", hjust = 0, size = 2.8, colour = "#D55E00") +
  annotate("text", x = 54, y = 82, label = "45 variants (7.8%)\n>50% reduction", hjust = 0, size = 2.8, colour = "#D55E00") +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                     sec.axis = dup_axis(breaks = top_breaks, labels = c("0", "-0.5", "-1", "-2", "-3.32"),
                                         name = expression(log[2](aFC)~cutoff))) +
  scale_y_continuous(limits = c(0, 285), breaks = seq(0, 250, 50), expand = expansion(mult = 0)) +
  labs(x = "Expression reduction threshold (%)", y = "Rare causal eQTL variants") + theme_classic(base_size = 9)
print(figure)

# %% Save aggregate outputs
#| Save the figure and threshold counts. There is no new significance test. The 105/574 value is at log2(aFC) < −0.5, or approximately 29.3% reduction. The >50% value is 45/574.
save_plot(figure, "supplements/S28", 130, 100)
readr::write_tsv(landmarks, output_file("supplements/S28_threshold_counts.tsv"))
