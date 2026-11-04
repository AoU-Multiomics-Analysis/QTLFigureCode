# Select the original rare all-assay rows and the comparator-gene rows.
# Run from the repository root. No input is written or downloaded.
source("preparation/common.R")
source("preparation/Fig5/functions.R")
config <- prep_init()
reference <- arrow::read_parquet(prep_input("input",
  "part-00000-006c7897-e219-4036-abbf-af4482a5177e-c000.snappy.parquet")) %>%
  as_tibble() %>% fig5_reference()
files <- c(eQTL = "COMB.expression_SusieMerged.annotated.ancestry_skew.tsv.gz",
  pQTL = "COMB.proteomics_SusieMerged.annotated.ancestry_skew.tsv.gz",
  sQTL = "COMB.splicing_SusieMerged..annotated.ancestry_skew.tsv.gz")
effect_symbols <- c("HEXA", "PMM2", "GALM", "GAA", "SMPD1", "CTH", "PLG", "CC2D1A", "ACY1")
effect_genes <- reference %>% filter(targetFromSourceId %in% effect_symbols) %>% pull(targetId)
audit <- list()
raw <- imap_dfr(files, function(file, assay) {
  message("Read all-assay input: ", assay)
  path <- prep_input("susie", file)
  columns <- c("variant", "molecular_trait_id", "gvs_max_af", "pip", "posterior_mean",
    "clinvar_classification", "cs_id", "consequence", "gvs_max_subpop", fig5_af_columns)
  data <- as_tibble(data.table::fread(path, select = columns, na.strings = "NA")) %>%
    mutate(gene = str_extract(molecular_trait_id, "ENSG[0-9]+"), assay = assay)
  rare <- data %>% filter(gvs_max_af < 0.01)
  audit[[assay]] <<- rare %>% summarise(assay = first(assay), rows = n(),
    variants = n_distinct(variant), missing_gene = sum(is.na(gene)),
    missing_cs = sum(is.na(cs_id)))
  stopifnot(!anyNA(rare$gene), !anyNA(rare$cs_id))
  data %>% filter(gvs_max_af < 0.01 | gene %in% effect_genes)
})
input_audit <- bind_rows(audit) %>%
  mutate(assay = recode(assay, eQTL = "Expression", pQTL = "Proteomics", sQTL = "Splicing"))
pairs <- fig5_select_pairs(raw, reference) %>%
  mutate(assay = recode(assay, eQTL = "Expression", pQTL = "Proteomics", sQTL = "Splicing")) %>%
  select(variant, gene, assay, pip, gvs_max_af, clinvar_classification) %>% distinct()
base <- "output/clingen_inheritance/rare_variant_count_bars_all_assays_unfiltered"
prep_write_tsv(input_audit, "paper", paste0(base, "/input_audit.tsv"))
prep_write_tsv(pairs, "paper", paste0(base, "/selected_pairs.tsv"))
prep_save_rds(list(raw = raw, reference = reference, pairs = pairs),
  "source", "Fig5/preparation/all_assay_context.rds")
message("Selected all-assay rows. No PIP or effect cutoff was added to A/B.")
