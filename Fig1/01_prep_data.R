# =============================================================================
# Figure 1 — data prep. Reads raw inputs ONCE, derives every per-panel data
# frame, and writes them to fig1/data/<panel>.rds for the panel scripts.
# Run from project root:  Rscript fig1/01_prep_data.R
#
# Expected input_data/ layout (see README "Data & inputs"):
#   reference/    gencode gtf, ancestry preds, COMB sample list
#   permutations/{eqtl,sqtl,pqtl}/   tensorQTL cis_qtl.txt.gz per ancestry
#   finemapping/{eqtl,sqtl,pqtl}/    SuSiE-merged tables
#   annotated/    COMB annotated SuSiE tables (for enrichment)
#   external/     GTEx / INTERVAL / TOPMed reference tables
# =============================================================================
source("fig1/00_setup.R")

saveit <- function(x, name) {
  saveRDS(x, file.path(DATA_DIR, paste0(name, ".rds")))
  message("  wrote ", file.path(DATA_DIR, paste0(name, ".rds")))
}

# Ancestry code from a fine-mapping filename (AoU_EAS_9k_... -> EAS, COMB kept).
clean_group <- function(x) str_remove_all(basename(x), "AoU_|_9k.*")

# =============================================================================
# Shared reference data (read once, used by several panels)
# =============================================================================
message("Loading gene annotations ...")
gene_annotations <- import(in_path("reference", "gencode.v48.annotation.gtf.gz"))
gene_types <- gene_annotations %>% data.frame() %>%
  select(gene_id, gene_type) %>% distinct() %>%
  mutate(gene_id = str_remove(gene_id, "\\..*"))

message("Loading ancestry counts ...")
RNA_COMB_sample_df <- fread(in_path("reference", "COMB_eQTL_9k_sample_list.tsv"))
ancestry_counts <- LoadAncestryCounts(
  in_path("reference", "echo_v4_r2.ancestry_preds.tsv"),
  RNA_COMB_sample_df$research_id
)

# =============================================================================
# QTL permutation results (Panels C, D)
# =============================================================================
# eQTL and sQTL permutation files share basenames ({ANCESTRY}.cis_qtl.txt.gz),
# so they are kept in separate subdirectories to avoid collisions.
read_permutations <- function(subdir, assay) {
  files <- list.files(in_path("permutations", subdir),
                      pattern = "\\.cis_qtl\\.txt\\.gz$", full.names = TRUE)
  purrr::map_dfr(files, function(f) {
    group <- basename(f) %>%
      str_remove("\\.cis_qtl\\.txt\\.gz$") %>%
      str_remove_all("AoU_10750_|_pQTL")
    fread(f) %>% mutate(group = group,
                        FDR = p.adjust(pval_beta, method = "fdr"),
                        Assay = assay)
  })
}

message("Loading QTL permutations ...")
QTL_data <- bind_rows(
  read_permutations("eqtl", "Expression"),
  read_permutations("pqtl", "Protein"),
  read_permutations("sqtl", "Splicing")
)

# =============================================================================
# Fine-mapping results (Panels E, F, G)
# =============================================================================
read_finemapping <- function(subdir, assay) {
  files <- list.files(in_path("finemapping", subdir),
                      pattern = "SusieMerged", full.names = TRUE)
  purrr::map_dfr(files, function(f) fread(f) %>% mutate(group = clean_group(f), Assay = assay))
}

message("Loading fine-mapping results ...")
eQTL_finemapping <- read_finemapping("eqtl", "Expression")
sQTL_finemapping <- read_finemapping("sqtl", "Splicing")
pQTL_susie       <- fread(in_path("finemapping", "pqtl", "pQTL_fm_data_10k.tsv"))

# =============================================================================
# Panel B — ancestry counts in cohort
# =============================================================================
saveit(ancestry_counts, "B")

# =============================================================================
# Panel C — proportion of genes with significant QTLs
# =============================================================================
combined_count <- sum(ancestry_counts$n)
prop_sig_data <- QTL_data %>%
  SummarizeSigQTLs(ancestry_counts) %>%
  left_join(AncestryColorPalette, by = c("group" = "Group")) %>%
  mutate(n = case_when(group == "AoU" ~ combined_count, TRUE ~ n))

add_aou_label <- function(df) {
  df %>% mutate(label = case_when(group == "AoU" ~ paste0("AoU", "\n(", num_sig, ")"), TRUE ~ group))
}
expression_only     <- prop_sig_data %>% filter(Assay == "Expression") %>% add_aou_label()
splicing_proteomics <- prop_sig_data %>% filter(Assay != "Expression") %>% add_aou_label()
saveit(list(expression_only = expression_only,
            splicing_proteomics = splicing_proteomics), "C")

# =============================================================================
# Panel D — novel eGene + sGene counts vs external consortia
# (pGene deferred: needs UKBB pQTL reference — see NOTES.md)
# =============================================================================
message("Loading external QTL references ...")
ExternaleQTLs <- LoadExternaleQTLs(
  in_path("external", "Whole_Blood.v10.eGenes.txt.gz"),
  in_path("external", "INTERVAL_eQTL_summary.xlsx"),
  in_path("external", "Whole_blood.maf0001.100.permutations.txt")
) %>% filter(Significant == TRUE)

ExternalsQTLs <- LoadExternalsQTLs(
  in_path("external", "Whole_Blood.v10.sGenes.txt.gz"),
  in_path("external", "INTERVAL_cis_sQTL_summary.xlsx")
) %>% filter(Significant == TRUE)

novel_eQTLs <- QTL_data %>%
  filter(group == "COMB" & Assay == "Expression") %>%
  mutate(phenotype_id = str_remove(phenotype_id, "\\..*")) %>%
  left_join(gene_types, by = c("phenotype_id" = "gene_id")) %>%
  filter(qval < .05) %>%
  ExtractNovelQTLs(ExternaleQTLs$gene_id) %>%
  mutate(type = "eGene")

novel_sQTLs <- QTL_data %>%
  filter(group == "COMB" & Assay == "Splicing") %>%
  mutate(phenotype_id = str_remove_all(phenotype_id, ".*:|\\..*")) %>%
  left_join(gene_types, by = c("phenotype_id" = "gene_id")) %>%
  filter(qval < .05) %>%
  ExtractNovelQTLs(ExternalsQTLs$gene_id) %>%
  mutate(type = "sGene")

gene_type_levels <- c("Protein-Coding", "lncRNA", "non-coding RNA")
novel_qtl_count <- bind_rows(novel_eQTLs, novel_sQTLs) %>%
  mutate(group = factor(group, levels = gene_type_levels),
         type  = factor(type,  levels = c("eGene", "sGene"))) %>%
  tidyr::complete(group, type, fill = list(count = 0)) %>%
  mutate(group = factor(group, levels = gene_type_levels),
         type  = factor(type,  levels = c("eGene", "sGene")))
saveit(novel_qtl_count, "D")

# =============================================================================
# Panel E — CDF of credible sets per eGene (COMB)
# =============================================================================
eQTL_cs_number <- eQTL_finemapping %>%
  mutate(Assay = "Expression") %>%
  distinct(molecular_trait_id, cs_id, group, Assay)
ccdf <- eQTL_cs_number %>% filter(Assay == "Expression") %>% ComputeCredibleSetCDF()
saveit(ccdf, "E")

# =============================================================================
# Panel F — credible-set sizes across pQTL ancestry groups
# =============================================================================
ref_means <- data.frame(group = "EUR", mean_cs = 9.7)   # UKBB pQTL reference average
pQTL_cs_sizes <- pQTL_susie %>%
  mutate(Assay = "Protein") %>%
  distinct(cs_id, cs_size, group, Assay)
common_var_means <- pQTL_susie %>%
  filter(!ALT_FREQS < .01 & !cs_size == 1) %>%
  mutate(Assay = "Protein") %>%
  distinct(cs_id, cs_size, group, Assay) %>%
  group_by(group) %>%
  summarize(mean_cs = median(cs_size)) %>%
  mutate(group = case_when(group == "COMB" ~ "AoU", TRUE ~ group))
saveit(list(pQTL_cs_sizes = pQTL_cs_sizes,
            common_var_means = common_var_means,
            ref_means = ref_means), "F")

# =============================================================================
# Panel G — proportion of credible sets with a single causal variant
# =============================================================================
finemapping_cs <- bind_rows(
  eQTL_finemapping %>% mutate(Assay = "Expression") %>% select(group, pip, Assay, cs_id),
  pQTL_susie       %>% mutate(Assay = "Protein")    %>% select(group, pip, Assay, cs_id),
  sQTL_finemapping %>% mutate(Assay = "Splicing")   %>% select(group, pip, Assay, cs_id)
)
summary_data <- finemapping_cs %>%
  distinct(cs_id, group, Assay) %>%
  dplyr::count(Assay, group) %>%
  dplyr::rename("num_cs" = "n")
causal_variant_count <- finemapping_cs %>%
  filter(pip > 0.9) %>%
  distinct(cs_id, group, Assay) %>%
  dplyr::count(Assay, group) %>%
  dplyr::rename("number_causal" = "n")
saveit(list(summary_data = summary_data,
            causal_variant_count = causal_variant_count), "G")

# =============================================================================
# Panel H — fine-mapping functional enrichment by assay
# =============================================================================
variant_functional_groups <- tibble::tribble(
  ~name, ~functional_group,
  # UTR
  "3_prime_utr_variant", "Regulatory Annotation",
  "5_prime_utr_variant", "Regulatory Annotation",
  # Regulatory elements
  "CTCF.bound", "Regulatory Annotation",
  "DNase.H3K4me3", "Regulatory Annotation",
  "FANTOM5", "Regulatory Annotation",
  "PLS", "Regulatory Annotation",
  "pELS", "Regulatory Annotation",
  "dELS", "Regulatory Annotation",
  # Coding
  "coding_sequence_variant", "Coding Annotation",
  "frameshift_variant", "Coding Annotation",
  "inframe_deletion", "Coding Annotation",
  "inframe_insertion", "Coding Annotation",
  "missense_variant", "Coding Annotation",
  "synonymous_variant", "Coding Annotation",
  "start_lost", "Coding Annotation",
  "stop_gained", "Coding Annotation",
  "stop_lost", "Coding Annotation",
  "stop_retained_variant", "Coding Annotation",
  # Splice
  "splice_acceptor_variant", "Splice Annotation",
  "splice_donor_variant", "Splice Annotation",
  "splice_donor_5th_base_variant", "Splice Annotation",
  "splice_donor_region_variant", "Splice Annotation",
  "splice_region_variant", "Splice Annotation",
  "splice_polypyrimidine_tract_variant", "Splice Annotation",
  # Noncoding RNA
  "non_coding_transcript_exon_variant", "Noncoding RNA",
  "non_coding_transcript_variant", "Noncoding RNA",
  "mature_mirna_variant", "Noncoding RNA",
  # Gene proximity / background
  "downstream_gene_variant", "Intronic/Intergenic",
  "upstream_gene_variant", "Intronic/Intergenic",
  "intron_variant", "Intronic/Intergenic",
  "intergenic_variant", "Intronic/Intergenic"
)

message("Computing fine-mapping annotation enrichment ...")
eQTL_variant_enrichment <- fread(in_path("annotated", "AoU_COMB_9k_eQTL_SusieMerged.annotated.tsv")) %>%
  ComputeFinemappingAnnotationEnrichment() %>% mutate(type = "eQTL")
sQTL_variant_enrichment <- fread(in_path("annotated", "AoU_COMB_9k_sQTL_SusieMerged.annotated.tsv")) %>%
  ComputeFinemappingAnnotationEnrichment() %>% mutate(type = "sQTL")
pQTL_variant_enrichment <- fread(in_path("annotated", "pQTL_susie_annotated.tsv")) %>%
  ComputeFinemappingAnnotationEnrichment() %>% mutate(type = "pQTL")
QTLFinemappingEnrichments <- bind_rows(
  pQTL_variant_enrichment, eQTL_variant_enrichment, sQTL_variant_enrichment
)
saveit(list(enrichment = QTLFinemappingEnrichments,
            functional_groups = variant_functional_groups), "H")

message("Data prep complete.")
