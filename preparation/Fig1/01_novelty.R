# Recompute current target-gene novelty from COMB tensorQTL and external results.
# Run from the repository root. All writes go to QTL_PREP_ROOT.
source("preparation/common.R")
source("preparation/Fig1/functions.R")
suppressPackageStartupMessages(library(readxl))
config <- prep_init()
base <- "output/qtl_circos_simplified_options"
current <- file.path(base, "N_main_figure_reconciled_20261027")
corrected <- file.path(base, "M_corrected_targets")
read_paper <- function(path, ...) read_tsv(prep_input("paper", path), show_col_types = FALSE, progress = FALSE, ...)
read_reference <- function(path, ...) read_tsv(prep_reference("paper", path), show_col_types = FALSE, progress = FALSE, ...)
read_input <- function(path, ...) read_tsv(prep_input("input", path), show_col_types = FALSE, progress = FALSE, ...)
write_current <- function(data, name) {
  prep_write_tsv(data, "paper", file.path(current, name))
  prep_write_tsv(data, "paper", file.path(corrected, name))
}

message("Figure 1: Read the gene reference prepared by fig1_gene_reference.")
genes <- read_paper(file.path(base, "gene_reference.tsv"))
hgnc <- read_paper("output/figure1_novel_gene_categories/annotations/hgnc_complete_set_2026-10-15.txt",
                  col_types = cols(.default = col_character()))
ukbb_file <- prep_input("input", "Fig1Data/41586_2023_6592_MOESM3_ESM.xlsx")
panel <- read_xlsx(ukbb_file, sheet = "ST3", skip = 2) %>%
  transmute(assay_target = `Assay Target`, gene_name = `Gene symbol`, uniprot_id = UniProt) %>%
  separate_longer_delim(gene_name, delim = ";") %>% mutate(gene_name = str_trim(gene_name))
st9 <- read_xlsx(ukbb_file, sheet = "ST9", skip = 4)
mapping <- fig1_ukbb_target_map(genes, hgnc, panel)
needed <- bind_rows(panel %>% select(gene_name), st9 %>% filter(`cis/trans` == "cis") %>%
  transmute(gene_name = `cis gene`) %>% separate_longer_delim(gene_name, delim = ",") %>%
  mutate(gene_name = str_trim(gene_name))) %>% distinct()
target_map <- mapping$map %>% semi_join(needed, by = "gene_name")
protein <- fig1_ukbb_cis_targets(st9, target_map)
measured <- panel %>% left_join(target_map, by = "gene_name", relationship = "many-to-many")
prep_write_tsv(target_map, "paper", file.path(corrected, "target_gene_map.tsv"))
prep_write_tsv(mapping$protein_identity, "paper", file.path(corrected, "uniprot_identity_mapping.tsv"))
prep_write_tsv(mapping$ambiguous %>% semi_join(needed, by = "gene_name"), "paper",
               file.path(corrected, "ambiguous_target_symbols.tsv"))
prep_write_tsv(protein$unresolved, "paper", file.path(corrected, "unmapped_ukbb_targets.tsv"))
prep_write_tsv(protein$audit, "paper", file.path(corrected, "ukbb_target_mapping_audit.tsv"))
prep_write_tsv(measured, "paper", file.path(corrected, "ukbb_panel_coverage.tsv"))

message("Figure 1: Read significant same-assay reference gene sets.")
gtex <- function(path, assay_name) read_input(file.path("Fig1Data", path),
  col_types = cols(.default = col_character())) %>% filter(as.double(qval) < .05) %>%
  transmute(assay = assay_name, gene_id = str_extract(gene_id, "ENSG[0-9]+"), study = "GTEx")
known <- bind_rows(
  gtex("Whole_Blood.v10.eGenes.txt.gz", "eQTL"),
  gtex("Whole_Blood.v10.sGenes.txt.gz", "sQTL"),
  read_xlsx(prep_input("input", "Fig1Data/INTERVAL_eQTL_summary.xlsx"), skip = 3) %>%
    filter(cis_eQTL_significant == 1) %>% transmute(assay = "eQTL", gene_id = feature_id, study = "INTERVAL"),
  read_xlsx(prep_input("input", "Fig1Data/INTERVAL_cis_sQTL_summary.xlsx"), skip = 3) %>%
    transmute(assay = "sQTL", gene_id = str_extract(gene_id, "ENSG[0-9]+"), study = "INTERVAL"),
  read_input("Fig1Data/Whole_blood.maf0001.100.permutations.txt",
             col_types = cols(.default = col_character())) %>% filter(as.double(qvalue) < .05) %>%
    transmute(assay = "eQTL", gene_id = str_extract(phenotype_id, "ENSG[0-9]+"), study = "Topmed"),
  protein$mapped %>% transmute(assay = "pQTL", gene_id, study = "UKBB-PPP")
) %>% filter(!is.na(gene_id), str_detect(gene_id, "^ENSG[0-9]+$")) %>% distinct()
known_compared <- prep_has_reference("paper", file.path(current, "paper_known_genes_by_study.tsv"))
if (known_compared) {
  saved_known <- read_reference(file.path(current, "paper_known_genes_by_study.tsv"))
  stopifnot(nrow(anti_join(known, saved_known, by = names(known))) == 0L,
            nrow(anti_join(saved_known, known, by = names(known))) == 0L)
}
write_current(known, "paper_known_genes_by_study.tsv")

message("Figure 1: Apply the source q-value cutoff and the current gene filters.")
raw_files <- c(eQTL = "eQTL/COMB.cis_qtl.txt.gz", sQTL = "sQTL/COMB.cis_qtl.txt.gz",
               pQTL = "pQTL/AoU_10750_COMB_pQTL.cis_qtl.txt.gz")
raw <- imap_dfr(raw_files, function(path, assay_name) {
  read_input(file.path("Fig1Data", path), col_types = cols(.default = col_character())) %>%
    transmute(assay = assay_name, phenotype_id, variant_id,
      source_pval_nominal = pval_nominal, source_pval_beta = pval_beta, source_qval = qval)
})
novelty <- fig1_classify_novelty(raw, genes, known, measured$gene_id)
stopifnot(nrow(novelty %>% filter(assay == "pQTL", gene_name %in% protein$unresolved$reference_symbol)) == 0L)
flags_compared <- prep_has_reference("paper", file.path(current, "paper_novelty_by_trait.tsv"))
stopifnot(nrow(novelty) == 44421L)
if (flags_compared) {
  saved <- read_reference(file.path(current, "paper_novelty_by_trait.tsv"),
                     col_types = cols(source_qval = col_character()))
comparison <- novelty %>% select(assay, phenotype_id, source_qval, gene_type, novel_for_figure) %>%
  full_join(saved %>% select(assay, phenotype_id, source_qval, gene_type, novel_for_figure),
    by = c("assay", "phenotype_id"), suffix = c(".rebuilt", ".saved"), relationship = "one-to-one")
stopifnot(nrow(comparison) == 44421L, !anyNA(comparison$novel_for_figure.rebuilt),
  all(comparison$novel_for_figure.rebuilt == comparison$novel_for_figure.saved),
  all(as.double(comparison$source_qval.rebuilt) == as.double(comparison$source_qval.saved)),
  all(coalesce(comparison$gene_type.rebuilt, "UNANNOTATED") == coalesce(comparison$gene_type.saved, "UNANNOTATED")))
}
counts <- novelty %>% filter(novel_for_figure) %>% count(group, type, name = "count") %>%
  complete(group = c("Protein-Coding", "lncRNA", "non-coding RNA"),
           type = c("eGene", "sGene", "pGene"), fill = list(count = 0L))
totals <- novelty %>% filter(novel_for_figure) %>% count(assay, name = "novel_genes")
stopifnot(identical(totals$novel_genes, c(10148L, 1143L, 2450L)),
          !anyDuplicated(novelty %>% filter(novel_for_figure) %>% select(assay, gene_id)))
write_current(novelty, "paper_novelty_by_trait.tsv")
write_current(counts, "novelty_counts_in_C_format.tsv")
write_current(novelty %>% filter(novel_for_figure) %>% count(assay, gene_type), "novel_genes_by_biotype.tsv")
coverage <- novelty %>% filter(assay == "pQTL", novel_for_figure) %>% count(novelty_context)
stopifnot(setequal(coverage$n, c(347L, 796L)))
write_current(coverage, "novel_pqtl_reference_coverage.tsv")
jsonlite::write_json(list(source_records = nrow(raw),
  saved_novelty_comparison = if (flags_compared) "passed" else "not_run",
  saved_reference_set_comparison = if (known_compared) "passed" else "not_run",
  all_44421_saved_calls_match = if (flags_compared) TRUE else NULL,
  all_source_q_values_match = if (flags_compared) TRUE else NULL,
  reference_set_differences = if (known_compared) 0L else NULL, totals = totals,
  pQTL_reference_coverage = coverage, source_q_values_recalculated = FALSE,
  rule = "Source qval < 0.05; same-assay reference target genes; current GENCODE biotype filters",
  caveat = "Target-gene novelty relative to the specified studies, not novel independent signals"),
  prep_output("paper", file.path(current, "novelty_provenance.json")), pretty = TRUE, auto_unbox = TRUE, null = "null")
message("Prepared Figure 1 novelty: 10,148 eGenes, 2,450 sGenes, 1,143 pGenes. Saved-call comparison: ",
        if (flags_compared) "passed" else "not run (reference absent)", ".")
