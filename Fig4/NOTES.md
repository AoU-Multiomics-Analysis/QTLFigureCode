# Figure 4 — notes and deferred work

## Result-affecting change: the duplicated conditional-TWAS table

`Fig4.R` computed `ColocAdjustedConditionalTWAS` **twice**. The first copy
(`Fig4.R:590`) referenced `fdr_rare_sig`, which is never defined
(`t_fdrboth = fdr_rare_cond < .1 & fdr_rare_sig`) — a typo for `rare_zscore_fdr`
— and used `p.adjust(method = "fdr")`. `GeneTWASData`, which Panel D depends on,
was built from that first copy. The second copy (`Fig4.R:636`) recomputed the
table cleanly with BH adjustment and a correct
`t_fdrboth = fdr_rare_cond < 0.10 & rare_zscore_fdr < 0.10`, but nothing used it.

`gene_twas_data()` keeps **only the clean (second) version** and builds
`GeneTWASData` from it. Panel D may therefore differ from a figure produced by
running the original top to bottom.

## GeneBayes was read but never used

`Fig4.R:91` read a GeneBayes s_het table from
`41588_2024_1820_MOESM4_ESM.xlsx` (sheet "Supplementary Table 1"), but nothing
downstream used it — every `shet` / `log_shet` comes from the `shet` column
already in `ConditionalTWAS`. The read is dropped, which also removes the
`readxl` dependency from Fig4. (This is the same supplement Fig1 (deferred) and
Fig3 reference.)

## Binary-predictor enrichment is computed but not plotted

Panel D fits enrichment for two predictor sets — `continuous`
(`max_rare`, `max_common`, `log_shet`, `log2(TPM)`) and `binary`
(`ever_rare`, `ever_common`, …) — but plots only the continuous `max_rare` /
`max_common` terms. Both are computed (`ENRICH_PREDICTORS` in `00_setup.R`); the
binary fits are available for inspection.

## Hand-curated heme/immune keyword list

The heme/immune gene set is defined by a regex of disease-name keywords
(`HEME_IMMUNE_KW` in `00_setup.R`) matched against ClinGen `diseaseFromSource`.
It is a curated list, not an ontology mapping.

## Inputs still missing a `gs://` origin

`inputs.tsv` records these local-only. The two `high_level_coloc_outputs/`
tables are extracted from `gs://fc-secure-…/high_level_coloc_outputs.tar.gz`
(README documents the staging command); `R2_TWAS_combined.tsv`,
`ConditionalTWAS.cleaned.tsv.gz`, `reference/interactions.tsv` and the ClinGen
parquet need URIs or pre-staging.
