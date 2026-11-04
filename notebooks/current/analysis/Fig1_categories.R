# %% Build the curated coding-gene categories
#| Filter the novel genes to protein-coding targets. Map curated annotation IDs before counting each category.
suppressPackageStartupMessages({library(arrow); library(jsonlite)})
gene_reference_all <- genes
genes <- novelty %>% filter(novel_for_figure, gene_type == "protein_coding") %>%
  distinct(assay, gene_id, gene_name)
parent <- "output/figure1_novel_gene_categories"
ann <- file.path(parent, "annotations")
out <- file.path(parent, "reconciled_20261027")
paper_path <- function(...) prep_input("paper", file.path(...))
totals <- genes %>% count(assay, name = "total_novel_coding")
reference <- gene_reference_all %>% distinct(gene_id, gene_name)

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

message("Mapped approved HGNC symbols and UniProt accessions.")

# %% Map positive human RNA-processing annotations
#| Exclude NOT and ND annotations. Prefer the UniProt mapping and retain stable gene IDs.
message("Map positive human RNA-processing annotations to stable gene IDs.")

gaf_names <- c("database", "product_id", "symbol", "qualifier", "go_id", "reference",
               "evidence", "with_from", "aspect", "name", "synonyms", "product_type",
               "taxon", "date", "assigned_by", "extension", "isoform")

gaf_path <- paper_path(out, "annotations/human_RNA_processing.gaf")

gaf <- read_tsv(gaf_path, comment = "!", col_names = gaf_names,
                col_types = cols(.default = col_character()), quote = "")


stopifnot(all(gaf$taxon == "taxon:9606"))

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

rna_genes <- rna %>% filter(!is.na(gene_id)) %>% distinct(gene_id)

message("Mapped positive RNA-processing annotations.")

# %% Map cancer, chromatin, disease and constraint annotations
#| Resolve cancer drivers and EpiFactors with HGNC IDs and unique GENCODE symbols.
#| Retain Moderate, Strong, or Definitive ClinGen evidence. Flagged or missing constraint values remain unknown.
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

tf <- read_lines(paper_path(ann, "human_tfs_ensembl_v1.01.txt")) %>% str_trim() %>%
  keep(~ str_detect(.x, "^ENSG[0-9]+$")) %>% unique()

clingen <- read_parquet(paper_path("output/clingen_inheritance/evidence_clingen.parquet")) %>%
  as_tibble() %>% filter(confidence %in% c("Moderate", "Strong", "Definitive")) %>% distinct(targetId)

constraint <- read_tsv(paper_path(ann, "gnomad_v2.1.1_by_gene.txt.bgz"),
  col_types = cols_only(gene_id = col_character(), oe_lof_upper = col_double(), constraint_flag = col_character()))

stopifnot(!anyDuplicated(constraint$gene_id), n_distinct(cancer$gene_id) == 591L)

message("Mapped curated gene annotations and loss-of-function constraint.")

# %% Join clinical drug stages to target mechanisms
#| Expand each drug mechanism to its target genes. Keep clinical or approved drugs with max_phase >= 1 and no withdrawal flag.
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

target_status <- mechanisms %>% group_by(gene_id) %>%
  summarise(approved_drug_target = any(approved_drug),
            clinical_or_approved_target = any(clinical_or_approved_drug),
            single_protein_approved_target = any(approved_drug & targetType == "single protein"),
            withdrawn_drug_target = any(withdrawn_drug),
            .groups = "drop")

drug <- target_status

stopifnot(!anyDuplicated(drug$gene_id))

message("Joined clinical drug stages to target mechanisms.")

# %% Join base category memberships to novel coding genes
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

message("Constructed the six base category memberships.")

# %% Add histone families and the matched GO annotation unions
#| Add HGNC histone families to the transcription/chromatin union.
#| For plasma membrane and signal transduction, include is_a and part_of descendants.
#| Keep positive human annotations. Exclude NOT and ND. Deduplicate stable gene IDs across roots.
local_genes <- genes %>% distinct(gene_id, gene_name)
id_counts <- hgnc %>% filter(!is.na(ensembl_gene_id)) %>% count(ensembl_gene_id, name = "id_matches")
id_map <- hgnc %>% inner_join(id_counts %>% filter(id_matches == 1L), by = "ensembl_gene_id") %>%
  transmute(gene_id = ensembl_gene_id, exact_hgnc_id = hgnc_id)
symbol_ids <- hgnc %>% transmute(gene_name = symbol, symbol_hgnc_id = hgnc_id)
identity <- local_genes %>% left_join(id_counts, by = c("gene_id" = "ensembl_gene_id")) %>%
  left_join(id_map, by = "gene_id") %>% left_join(symbol_ids, by = "gene_name") %>%
  mutate(hgnc_id = case_when(id_matches == 1L ~ exact_hgnc_id,
                             is.na(id_matches) ~ symbol_hgnc_id, TRUE ~ NA_character_)) %>%
  left_join(hgnc %>% select(hgnc_id, gene_group), by = "hgnc_id", relationship = "many-to-one")
histones <- identity %>% separate_longer_delim(gene_group, delim = "|") %>%
  filter(str_detect(str_trim(gene_group), "^H[1-4](A|B)? histones$")) %>% distinct(gene_id)
hgnc_genes <- identity %>% filter(!is.na(hgnc_id)) %>% distinct(hgnc_id, gene_id)
protein_hgnc <- hgnc %>% separate_longer_delim(uniprot_ids, delim = "|") %>%
  filter(!is.na(uniprot_ids), uniprot_ids != "") %>% distinct(uniprot_id = uniprot_ids, hgnc_id)
unique_reference_symbols <- reference %>% group_by(gene_name) %>%
  filter(n_distinct(gene_id) == 1L) %>% ungroup() %>%
  transmute(symbol = gene_name, fallback_gene_id = gene_id)
message("Mapped HGNC histone families.")

# %% Read GO terms and construct descendant sets
ontology_path <- paper_path("output/uncategorized_novel_egenes_20261027/annotations/go-basic_GO_release_2026-08-05.json.gz")
gaf_path <- paper_path("output/uncategorized_novel_egenes_20261027/annotations/HUMAN-uniprot_GO_release_2026-08-05.gaf.gz")
graph <- jsonlite::fromJSON(gzfile(ontology_path), simplifyVector = FALSE)$graphs[[1]]
go_id <- function(x) str_replace(str_extract(x, "GO[_:][0-9]{7}$"), "_", ":")
go_nodes <- tibble(
  go_id = map_chr(graph$nodes, ~go_id(.x$id)),
  type = map_chr(graph$nodes, ~.x$type %||% ""),
  label = map_chr(graph$nodes, ~.x$lbl %||% ""),
  deprecated = map_lgl(graph$nodes, ~isTRUE(.x$meta$deprecated))) %>%
  filter(!is.na(go_id), type == "CLASS", !deprecated)
valid_go_nodes <- graph$nodes[map_chr(graph$nodes, ~go_id(.x$id)) %in% go_nodes$go_id]
go_alternatives <- map_dfr(valid_go_nodes, function(node) {
  properties <- node$meta$basicPropertyValues
  if (!length(properties)) return(tibble())
  map_dfr(properties, function(item) {
    if (!str_ends(item$pred, "#hasAlternativeId")) return(tibble())
    tibble(alternative = go_id(item$val), primary = go_id(node$id))
  })
})
go_edges <- tibble(
  child = map_chr(graph$edges, ~go_id(.x$sub)),
  parent = map_chr(graph$edges, ~go_id(.x$obj)),
  predicate = map_chr(graph$edges, ~.x$pred))
go_roots <- tibble(root = c("GO:0005886", "GO:0007165"),
                   category = c("plasma_membrane", "signal_transduction"))
stopifnot(str_detect(graph$meta$version, "/2026-07-26/"),
          go_nodes$label[match(go_roots$root, go_nodes$go_id)] == c("plasma membrane", "signal transduction"))
descendants <- go_roots %>% mutate(go_id = map(root, ~fig1_go_descendants(go_nodes, go_edges, .x))) %>%
  unnest(go_id)
message("Constructed plasma-membrane and signal-transduction GO descendant sets.")

# %% Map positive GO annotations and form gene unions
membrane_gaf <- read_tsv(gaf_path, comment = "!", col_names = gaf_names,
  col_types = cols(.default = col_character()), quote = "", progress = FALSE) %>%
  filter(database == "UniProtKB", str_detect(taxon, "^taxon:9606($|\\|)"),
         !str_detect(coalesce(qualifier, ""), "(^|\\|)NOT(\\||$)"), evidence != "ND") %>%
  left_join(go_alternatives, by = c("go_id" = "alternative"), relationship = "many-to-one") %>%
  mutate(go_id = coalesce(primary, go_id), uniprot_id = str_remove(product_id, "-[0-9]+$")) %>%
  inner_join(descendants, by = "go_id", relationship = "many-to-many")
# HGNC UniProt matches take precedence, including matches outside the local gene list.
# Use the approved symbol only when there is no protein-accession match.
membrane_evidence <- membrane_gaf %>%
  left_join(protein_hgnc %>% rename(protein_hgnc_id = hgnc_id), by = "uniprot_id", relationship = "many-to-many") %>%
  left_join(hgnc %>% transmute(symbol, approved_hgnc_id = hgnc_id), by = "symbol", relationship = "many-to-one") %>%
  mutate(hgnc_id = coalesce(protein_hgnc_id, approved_hgnc_id)) %>%
  left_join(hgnc_genes, by = "hgnc_id", relationship = "many-to-many") %>%
  left_join(unique_reference_symbols, by = "symbol", relationship = "many-to-one") %>%
  mutate(gene_id = if_else(is.na(hgnc_id), fallback_gene_id, gene_id)) %>%
  semi_join(local_genes, by = "gene_id")
membrane_union <- membrane_evidence %>% distinct(gene_id, category) %>%
  mutate(member = TRUE) %>% pivot_wider(names_from = category, values_from = member, values_fill = FALSE)
membership <- membership %>%
  mutate(transcription_chromatin_previous = transcription_chromatin,
         histone = gene_id %in% histones$gene_id,
         transcription_chromatin = transcription_chromatin | histone) %>%
  left_join(membrane_union, by = "gene_id", relationship = "many-to-one") %>%
  mutate(across(c(plasma_membrane, signal_transduction), ~replace_na(.x, FALSE)),
         plasma_membrane_signaling = plasma_membrane | signal_transduction)
message("Joined histone and GO union membership to novel coding genes.")

# %% Count the seven category memberships
categories <- c("cancer_driver", "mendelian_disease", "lof_constrained", "transcription_chromatin",
                "rna_processing_splicing", "plasma_membrane_signaling", "drug_target")
category_definitions <- definitions %>%
  bind_rows(tibble(category = "plasma_membrane_signaling", label = "Plasma membrane / signal transduction",
    definition = "Union of positive human GO:0005886 and GO:0007165 annotations and is_a/part_of descendants; exclude NOT and ND")) %>%
  mutate(definition = if_else(category == "transcription_chromatin",
    paste0(definition, "; HGNC H1, H2A, H2B, H3 and H4 histone families"), definition))
category_result <- fig1_category_counts(membership, categories)
category_counts <- category_result$counts %>% left_join(totals, by = "assay") %>%
  left_join(category_definitions %>% select(category, label), by = "category")
category_genes <- category_result$lists
stopifnot(nrow(category_counts) == 21L, n_distinct(local_genes$gene_id) == 4146L)
message("Calculated all seven category memberships and 21 assay/category counts.")
