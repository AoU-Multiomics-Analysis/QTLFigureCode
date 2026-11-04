# %% Assign ancestry groups and calculate the Figure 5 summaries
#| Count distinct variants within each inheritance group. Recalculate ClinVar coverage, the AR funnel, effect comparisons, and Fisher tests.
frequencies <- raw %>% select(variant, all_of(fig5_af_columns)) %>% distinct()
assignments <- fig5_variant_assignments(pairs, frequencies, reference)
summary <- fig5_counts_and_classes(assignments, pairs)
funnel <- fig5_funnel(raw, reference)
effects <- fig5_effect_panels(raw, funnel, reference)
coverage_tests <- fig5_coverage_tests(summary$coverage, assignments)
primary <- bind_rows(
  fig5_subset_test(funnel$sets$rare_all, funnel$sets$pip09,
    "pip09", "rare_all", "PIP > 0.9 versus starting pool"),
  fig5_subset_test(funnel$sets$rare_all, funnel$sets$posterior_lt_neg03,
    "posterior_lt_neg03", "rare_all", "PIP > 0.9 and effect < -0.3 versus starting pool")
)
incremental <- fig5_subset_test(funnel$sets$pip09, funnel$sets$posterior_lt_neg03,
  "posterior_lt_neg03", "pip09", "Additional effect filter within high-PIP variants")

stopifnot(n_distinct(assignments$variant) == 4479L, nrow(effects$D) == 14L, nrow(effects$E) == 16L)
message("Calculated Figure 5 counts, filters, effects, and statistical comparisons.")
