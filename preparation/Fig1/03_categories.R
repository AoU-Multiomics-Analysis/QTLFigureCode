#!/usr/bin/env Rscript
# Use one corrected target-gene set for the main figure and manuscript counts.
source("preparation/common.R")
source("preparation/Fig1/functions.R")
suppressPackageStartupMessages({library(arrow); library(jsonlite)})
config <- prep_init()
parent <- "output/figure1_novel_gene_categories"
current <- file.path(parent, "seven_categories_20261027")
paper_path <- function(...) prep_input("paper", file.path(...))
out <- "output/figure1_novel_gene_categories/reconciled_20261027"
old <- "output/figure1_novel_gene_categories"
ann <- file.path(old, "annotations")
novel_path <- paper_path("output/qtl_circos_simplified_options/M_corrected_targets/paper_novelty_by_trait.tsv")


message("Read corrected novel protein-coding targets and verify assay totals.")
genes <- read_tsv(novel_path, col_types = cols(.default = col_character())) %>%
  filter(novel_for_figure == "TRUE", gene_type == "protein_coding") %>%
  distinct(assay, gene_id, gene_name)
totals <- genes %>% count(assay, name = "total_novel_coding") %>%
  arrange(match(assay, c("eQTL", "sQTL", "pQTL")))
stopifnot(identical(totals$total_novel_coding, c(1207L, 2035L, 1143L)),
          !anyDuplicated(genes[c("assay", "gene_id")]))
write_tsv(totals, prep_output("paper", file.path(out, "assay_totals.tsv")))

reference <- read_tsv(paper_path("output/qtl_circos_simplified_options/gene_reference.tsv"),
                      col_types = cols(.default = col_character())) %>%
  distinct(gene_id, gene_name)
hgnc <- read_tsv(paper_path(ann, "hgnc_complete_set_2026-10-15.txt"),
                 col_types = cols(.default = col_character())) %>% filter(status == "Approved")
symbol_map <- reference %>% group_by(gene_name) %>%
  filter(n_distinct(gene_id) == 1L) %>% ungroup() %>%
  transmute(symbol = gene_name, gene_id, priority = 3L, method = "Unique GENCODE symbol")
approved <- hgnc %>% filter(!is.na(ensembl_gene_id)) %>%
  transmute(symbol, gene_id = ensembl_gene_id, priority = 2L, method = "HGNC approved symbol")
previous <- hgnc %>% filter(!is.na(ensembl_gene_id)) %>%
  separate_longer_delim(prev_symbol, delim = "|") %>% filter(!is.na(prev_symbol)) %>%
  transmute(symbol = prev_symbol, gene_id = ensembl_gene_id, priority = 4L, method = "HGNC previous symbol")
aliases <- hgnc %>% filter(!is.na(ensembl_gene_id)) %>%
  separate_longer_delim(alias_symbol, delim = "|") %>% filter(!is.na(alias_symbol)) %>%
  transmute(symbol = alias_symbol, gene_id = ensembl_gene_id, priority = 5L, method = "HGNC unique alias")
symbols <- bind_rows(approved, symbol_map, previous, aliases) %>% distinct() %>%
  group_by(symbol) %>% filter(priority == min(priority), n_distinct(gene_id[priority == min(priority)]) == 1L) %>%
  ungroup() %>% distinct(symbol, gene_id, method)
stopifnot(!anyDuplicated(symbols$symbol))
uniprot <- hgnc %>% filter(!is.na(ensembl_gene_id)) %>%
  separate_longer_delim(uniprot_ids, delim = "|") %>% filter(!is.na(uniprot_ids)) %>%
  group_by(uniprot_ids) %>% filter(n_distinct(ensembl_gene_id) == 1L) %>% ungroup() %>%
  distinct(uniprot_id = uniprot_ids, uniprot_gene_id = ensembl_gene_id)

message("Map positive human RNA-processing annotations to stable gene IDs.")
gaf_names <- c("database", "product_id", "symbol", "qualifier", "go_id", "reference",
               "evidence", "with_from", "aspect", "name", "synonyms", "product_type",
               "taxon", "date", "assigned_by", "extension", "isoform")
gaf_path <- paper_path(out, "annotations/human_RNA_processing.gaf")
gaf <- read_tsv(gaf_path, comment = "!", col_names = gaf_names,
                col_types = cols(.default = col_character()), quote = "")
expected_rows <- read_json(paper_path(out, "annotations/download_manifest.json"))$annotation_rows
stopifnot(nrow(gaf) == expected_rows, all(gaf$taxon == "taxon:9606"))
rna <- gaf %>%
  filter(database == "UniProtKB", !str_detect(coalesce(qualifier, ""), "(^|\\|)NOT(\\||$)"),
         evidence != "ND") %>%
  mutate(uniprot_id = str_remove(product_id, "-[0-9]+$")) %>%
  left_join(uniprot, by = "uniprot_id", relationship = "many-to-one") %>%
  left_join(symbols %>% rename(symbol_gene_id = gene_id, symbol_method = method),
            by = "symbol", relationship = "many-to-one") %>%
  mutate(gene_id = coalesce(uniprot_gene_id, symbol_gene_id),
         mapping = case_when(!is.na(uniprot_gene_id) ~ "HGNC UniProt ID",
                             !is.na(symbol_gene_id) ~ symbol_method, TRUE ~ "Unmapped"))
write_tsv(rna, prep_output("paper", file.path(out, "rna_processing_annotation_mapping.tsv")))
rna %>% distinct(product_id, symbol, gene_id, mapping) %>%
  write_tsv(prep_output("paper", file.path(out, "rna_processing_gene_mapping.tsv")))
rna_genes <- rna %>% filter(!is.na(gene_id)) %>% distinct(gene_id)

message("Reapply the original curated annotations to the corrected target list.")
curated_symbol_map <- reference %>% group_by(gene_name) %>%
  filter(n_distinct(gene_id) == 1L) %>% ungroup() %>%
  transmute(symbol = gene_name, reference_gene_id = gene_id)
entrez_map <- NULL
cancer_raw <- read_tsv(paper_path(ann, "ncg_v7.2_cancer_drivers.tsv"),
                       col_types = cols(.default = col_character())) %>%
  filter(type == "Canonical Cancer Driver") %>% distinct(entrez, symbol)
stopifnot(n_distinct(cancer_raw$entrez) == 591L)
entrez_map <- hgnc %>% filter(!is.na(entrez_id)) %>%
  group_by(entrez_id) %>% filter(n_distinct(ensembl_gene_id, na.rm = TRUE) == 1L) %>%
  ungroup() %>% distinct(entrez_id, ensembl_gene_id)
cancer <- cancer_raw %>%
  left_join(entrez_map, by = c("entrez" = "entrez_id"), relationship = "many-to-one") %>%
  left_join(curated_symbol_map, by = "symbol", relationship = "many-to-one") %>%
  mutate(gene_id = coalesce(ensembl_gene_id, reference_gene_id),
         mapping = case_when(!is.na(ensembl_gene_id) ~ "HGNC Entrez-to-Ensembl",
                             !is.na(reference_gene_id) ~ "Unique GENCODE symbol",
                             TRUE ~ "Unmapped"))

epi_raw <- read_csv(paper_path(ann, "epifactors_v2.1_proteins.csv"),
                    col_types = cols(.default = col_character()), na = c("", "#"))
epi <- epi_raw %>%
  mutate(hgnc_id = if_else(str_detect(HGNC_ID, "^HGNC:"), HGNC_ID, paste0("HGNC:", HGNC_ID))) %>%
  left_join(hgnc %>% select(hgnc_id, ensembl_gene_id), by = "hgnc_id", relationship = "many-to-one") %>%
  left_join(curated_symbol_map, by = c("HGNC_symbol" = "symbol"), relationship = "many-to-one") %>%
  mutate(gene_id = coalesce(ensembl_gene_id, reference_gene_id),
         mapping = case_when(!is.na(ensembl_gene_id) ~ "HGNC ID-to-Ensembl",
                             !is.na(reference_gene_id) ~ "Unique GENCODE symbol",
                             TRUE ~ "Unmapped"),
         # Include the named functions and their cofactors, with overlap allowed.
         chromatin_remodeling = str_detect(coalesce(Function, ""), "(^|, )Chromatin remodeling( cofactor)?(,|$)"),
         histone_dna_modification = str_detect(coalesce(Function, ""), "(^|, )(Histone modification|DNA modification)"))
write_tsv(cancer, prep_output("paper", file.path(out, "cancer_mapping_audit.tsv")))
write_tsv(epi, prep_output("paper", file.path(out, "epifactors_mapping_audit.tsv")))

tf <- read_lines(paper_path(ann, "human_tfs_ensembl_v1.01.txt")) %>% str_trim() %>%
  keep(~ str_detect(.x, "^ENSG[0-9]+$")) %>% unique()
clingen <- read_parquet(paper_path("output/clingen_inheritance/evidence_clingen.parquet")) %>%
  as_tibble() %>% filter(confidence %in% c("Moderate", "Strong", "Definitive")) %>% distinct(targetId)
write_tsv(clingen, prep_output("paper", file.path(out, "clingen_qualifying_gene_ids.tsv")))
constraint <- read_tsv(paper_path(ann, "gnomad_v2.1.1_by_gene.txt.bgz"),
  col_types = cols_only(gene_id = col_character(), oe_lof_upper = col_double(), constraint_flag = col_character()))
stopifnot(!anyDuplicated(constraint$gene_id), n_distinct(cancer$gene_id) == 591L)
drug_annotations <- file.path(parent, "drug_target_counts/annotations")
phase_path <- paper_path(drug_annotations, "chembl_clinical_molecules_2026-10-16.tsv")
mechanism_path <- paper_path(drug_annotations, "drug_mechanism_of_action_26.09.parquet")
molecule_path <- paper_path(drug_annotations, "drug_molecule_26.09.parquet")
phases <- read_tsv(phase_path, col_types = cols(.default = col_character())) %>%
  transmute(drug_id = molecule_chembl_id, chembl_name = pref_name,
            chembl_max_phase = as.double(max_phase),
            chembl_withdrawn = str_to_lower(withdrawn_flag) == "true")
drugs <- read_parquet(molecule_path) %>%
  select(drug_id = id, drug_name = name, maximumClinicalStage, drugType) %>%
  left_join(phases, by = "drug_id", relationship = "one-to-one")
stopifnot(!anyDuplicated(drugs$drug_id))
message("Expand curated drug mechanisms to distinct target gene and drug pairs.")
mechanisms <- read_parquet(mechanism_path) %>%
  select(mechanismOfAction, actionType, targetName, targetType, chemblIds, targets) %>%
  mutate(chemblIds = as.list(chemblIds), targets = as.list(targets)) %>%
  unnest_longer(chemblIds, values_to = "drug_id") %>%
  unnest_longer(targets, values_to = "gene_id") %>%
  filter(!is.na(gene_id), str_detect(gene_id, "^ENSG[0-9]+$")) %>%
  left_join(drugs, by = "drug_id", relationship = "many-to-one") %>% distinct()
stopifnot(!anyNA(mechanisms$maximumClinicalStage))
clinical_stages <- c("EARLY_PHASE_1", "PHASE_1", "PHASE_1_2", "PHASE_2", "PHASE_2_3",
                     "PHASE_3", "PREAPPROVAL", "PHASE_4", "APPROVAL")
stage_counts <- drugs %>% count(maximumClinicalStage, name = "drugs")

mechanisms <- mechanisms %>%
  mutate(opentargets_approval_stage = maximumClinicalStage == "APPROVAL",
         approved_drug = coalesce(chembl_max_phase == 4 & !chembl_withdrawn, FALSE),
         clinical_or_approved_drug = coalesce(chembl_max_phase >= 1 & !chembl_withdrawn, FALSE),
         withdrawn_drug = coalesce(chembl_withdrawn, FALSE))
write_tsv(mechanisms %>% filter(opentargets_approval_stage != approved_drug) %>%
            distinct(drug_id, drug_name, maximumClinicalStage, chembl_max_phase, chembl_withdrawn),
          prep_output("paper", file.path(out, "approval_annotation_discrepancies.tsv")))
target_status <- mechanisms %>% group_by(gene_id) %>%
  summarise(approved_drug_target = any(approved_drug),
            clinical_or_approved_target = any(clinical_or_approved_drug),
            single_protein_approved_target = any(approved_drug & targetType == "single protein"),
            withdrawn_drug_target = any(withdrawn_drug),
            .groups = "drop")
drug <- target_status
write_tsv(target_status, prep_output("paper", file.path(parent, "drug_target_counts/all_opentargets_drug_target_status.tsv")))
stopifnot(!anyDuplicated(drug$gene_id))
membership <- genes %>%
  left_join(constraint, by = "gene_id", relationship = "many-to-one") %>%
  left_join(drug, by = "gene_id", relationship = "many-to-one") %>%
  mutate(cancer_driver = gene_id %in% cancer$gene_id,
         mendelian_disease = gene_id %in% clingen$targetId,
         lof_constrained = if_else(!is.na(oe_lof_upper) & is.na(constraint_flag), oe_lof_upper < 0.35, NA),
         transcription_factor = gene_id %in% tf,
         chromatin_remodeling = gene_id %in% epi$gene_id[epi$chromatin_remodeling %in% TRUE],
         histone_dna_modification = gene_id %in% epi$gene_id[epi$histone_dna_modification %in% TRUE],
         transcription_chromatin = transcription_factor | chromatin_remodeling | histone_dna_modification,
         rna_processing_splicing = gene_id %in% rna_genes$gene_id,
         across(ends_with("target"), ~ replace_na(.x, FALSE)),
         drug_target = clinical_or_approved_target)
categories <- c("cancer_driver", "mendelian_disease", "lof_constrained", "transcription_chromatin",
                "rna_processing_splicing", "drug_target")
definitions <- tribble(
  ~category, ~label, ~definition,
  "cancer_driver", "Cancer drivers", "NCG 7.2 canonical cancer drivers",
  "mendelian_disease", "Mendelian disease genes", "ClinGen Moderate, Strong or Definitive gene-disease evidence; Open Targets 26.06",
  "lof_constrained", "LoF-constrained genes", "gnomAD v2.1.1 LOEUF <0.35 with no constraint flag; missing or flagged values remain unknown",
  "transcription_chromatin", "Transcription and chromatin regulation", "Union of Human TFs v1.01 and EpiFactors 2.1 chromatin remodeling and histone/DNA modification; stable gene ID deduplication",
  "rna_processing_splicing", "RNA processing and splicing", "Positive human GO:0006396 annotations, including is_a/part_of descendants and RNA splicing GO:0008380; all evidence except ND; exclude NOT; QuickGO snapshot 2026-10-27",
  "drug_target", "Drug targets", "Open Targets 26.09 drug mechanisms linked to ChEMBL 37 max_phase >=1 and withdrawn_flag=false; approved or clinical, including mapped protein family/complex members; snapshot 2026-10-16"
)
write_tsv(definitions, prep_output("paper", file.path(out, "category_definitions.tsv")))
write_tsv(membership, prep_output("paper", file.path(out, "novel_coding_gene_annotations.tsv")))
long <- membership %>% pivot_longer(all_of(categories), names_to = "category", values_to = "member")
counts <- long %>% group_by(assay, category) %>%
  summarise(n_genes = n_distinct(gene_id[member %in% TRUE]), unknown = sum(is.na(member)), .groups = "drop") %>%
  left_join(totals, by = "assay") %>% left_join(definitions, by = "category") %>%
  arrange(match(category, categories), match(assay, c("eQTL", "sQTL", "pQTL")))
lists <- long %>% filter(member %in% TRUE) %>% distinct(assay, category, gene_id, gene_name)
distinct_counts <- lists %>% group_by(category) %>% summarise(n_genes = n_distinct(gene_id), .groups = "drop") %>%
  left_join(definitions, by = "category") %>% arrange(match(category, categories))
write_tsv(counts, prep_output("paper", file.path(out, "category_counts_by_assay.tsv")))
write_tsv(lists, prep_output("paper", file.path(out, "category_gene_lists.tsv")))
write_tsv(distinct_counts, prep_output("paper", file.path(out, "distinct_category_counts.tsv")))
lists %>% distinct(category, gene_id, gene_name) %>% arrange(category, gene_id) %>%
  write_tsv(prep_output("paper", file.path(out, "distinct_category_gene_lists.tsv")))

# Expand the six-category membership with histones and the matched GO unions.
# The Python implementation is the archived scientific transform with portable roots.
raw_inputs <- c(file.path(parent, "annotations/hgnc_complete_set_2026-10-15.txt"),
  "output/uncategorized_novel_egenes_20261027/annotations/go-basic_GO_release_2026-08-05.json.gz",
  "output/uncategorized_novel_egenes_20261027/annotations/HUMAN-uniprot_GO_release_2026-08-05.gaf.gz")
walk(raw_inputs, ~prep_input("paper", .x))
expanded_files <- c("gene_identity_mapping.tsv", "histone_gene_mapping.tsv",
  "membrane_signaling_annotation_evidence.tsv.gz", "membrane_signaling_gene_union.tsv",
  "novel_coding_gene_annotations.tsv", "category_definitions.tsv", "category_counts_by_assay.tsv",
  "distinct_category_counts.tsv", "category_gene_lists.tsv", "distinct_category_gene_lists.tsv",
  "coverage_summary.tsv", "union_deduplication_audit.tsv", "assay_totals.tsv",
  "manuscript_sentence.txt", "data_audit.json")
walk(expanded_files, ~prep_output("paper", file.path(current, .x)))
status <- system2(Sys.getenv("PYTHON", "python3"), c("preparation/Fig1/expand_categories.py",
  "--paper-root", shQuote(config$paper_root), "--prepared-paper-root", shQuote(file.path(config$output_root, "paper"))))
if (status != 0L) stop("Figure 1 seven-category preparation failed.")
new_counts <- read_tsv(file.path(config$output_root, "paper", current, "category_counts_by_assay.tsv"), show_col_types = FALSE)
saved_counts <- read_tsv(prep_reference("paper", file.path(current, "category_counts_by_assay.tsv")), show_col_types = FALSE)
keys <- c("assay", "category", "n_genes", "total_novel_coding")
stopifnot(nrow(anti_join(new_counts, saved_counts, by = keys)) == 0L,
          nrow(anti_join(saved_counts, new_counts, by = keys)) == 0L)
new_lists <- read_tsv(file.path(config$output_root, "paper", current, "category_gene_lists.tsv"), show_col_types = FALSE)
saved_lists <- read_tsv(prep_reference("paper", file.path(current, "category_gene_lists.tsv")), show_col_types = FALSE)
stopifnot(nrow(anti_join(new_lists, saved_lists, by = names(saved_lists))) == 0L,
          nrow(anti_join(saved_lists, new_lists, by = names(saved_lists))) == 0L)
message("Figure 1 categories verified: all 21 count cells and gene memberships match the current seven-category data.")
