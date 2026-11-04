# %% Read GENCODE genes and correct UKBB protein identities
#| Stable gene IDs come from the GENCODE v48 GTF. HGNC symbols and UniProt accessions correct UKBB targets.
suppressPackageStartupMessages(library(readxl))
genes <- fig1_read_gene_gtf(prep_input("input", "gencode.v48.annotation.gtf.gz"))
stopifnot(nrow(genes) == 78686L)
read_paper <- function(path, ...) read_tsv(prep_input("paper", path), show_col_types = FALSE, progress = FALSE, ...)
read_input <- function(path, ...) read_tsv(prep_input("input", path), show_col_types = FALSE, progress = FALSE, ...)
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

message("Mapped UKBB protein targets with HGNC and UniProt identities.")

# %% Join significant blood-study targets to the source QTL results
#| Each reference uses its reported significant targets. Main-figure novelty uses source qval < 0.05.
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

stopifnot(nrow(novelty) == 44421L)
message("Classified all 44,421 source traits with assay-specific reference targets.")
