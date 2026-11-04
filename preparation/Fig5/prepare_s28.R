# Rebuild the historical 574-variant eQTL set and reduction thresholds.
source("preparation/common.R")
source("preparation/Fig5/functions.R")
config <- prep_init()
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
curve <- fig5_threshold_tables(selected, seq(0, 100, by = 1))
landmarks <- fig5_threshold_tables(selected,
  sort(unique(c(seq(0, 100, by = 5), 100 * (1 - 2^-0.5), 99))))
base <- "Supplemental Figures/01_Supplemental_drafts/Figure5_rare_eQTL_aFC_thresholds"
for (name in c("threshold_curve", "threshold_counts")) {
  observed <- if (name == "threshold_curve") curve$counts else landmarks$counts
  expected <- read_tsv(prep_reference("paper", paste0(base, "/", name, ".tsv")),
    show_col_types = FALSE)
  stopifnot(isTRUE(all.equal(observed, expected, check.attributes = FALSE, tolerance = 1e-10)))
  prep_write_tsv(observed, "paper", paste0(base, "/", name, ".tsv"))
}
prep_write_tsv(selected, "paper",
  "output/clingen_inheritance/reproduce_574_original_opentargets_selected.tsv")
prep_write_tsv(curve$variants, "paper",
  "output/clingen_inheritance/rare_eqtl_afc_thresholds/v2/distinct_variant_effects.tsv")
message("S28 matched: 574 variants; 45 exceed 50% reduction and 105 pass log2(aFC) < -0.5.")
