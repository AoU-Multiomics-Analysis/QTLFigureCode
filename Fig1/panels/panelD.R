# Panel D — novel eGene + sGene + pGene counts vs external consortia
# Inputs: perm (COMB expression/splicing/protein), gene_gtf, external references,
#         ukbb_pqtl
source("fig1/00_setup.R")

QTL_data <- load_input("perm")

gene_types <- import(input_path("gene_gtf")) %>% data.frame() %>%
  select(gene_id, gene_type) %>% distinct() %>%
  mutate(gene_id = str_remove(gene_id, "\\..*"))

ExternaleQTLs <- LoadExternaleQTLs(
  input_path("gtex_egenes"),
  input_path("interval_eqtl"),
  input_path("topmed_eqtl")
) %>% filter(Significant == TRUE)

ExternalsQTLs <- LoadExternalsQTLs(
  input_path("gtex_sgenes"),
  input_path("interval_sqtl")
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

# Novel pGenes vs the UKBB-PPP cis-pQTL reference. The original notebook
# referenced UKBB_pqtl_path but never assigned it, so this bar was silently
# absent; see NOTES.md. lncRNA is dropped because the proteomics assay cannot
# measure it, so a "novel lncRNA pGene" count would be an artefact of the join.
UKBB_pQTLs <- LoadUKBBpQTL(input_path("ukbb_pqtl"))

novel_pQTLs <- QTL_data %>%
  filter(group == "COMB" & Assay == "Protein") %>%
  filter(qval < .05) %>%
  mutate(phenotype_id = str_remove_all(phenotype_id, ".*_|\\..*")) %>%
  left_join(gene_types, by = c("phenotype_id" = "gene_id")) %>%
  ExtractNovelQTLs(UKBB_pQTLs$gene_id) %>%
  filter(group != "lncRNA") %>%
  mutate(type = "pGene")

gene_type_levels <- c("Protein-Coding", "lncRNA", "non-coding RNA")
type_levels      <- c("eGene", "sGene", "pGene")
novel_qtl_count <- bind_rows(novel_eQTLs, novel_sQTLs, novel_pQTLs) %>%
  mutate(group = factor(group, levels = gene_type_levels),
         type  = factor(type,  levels = type_levels)) %>%
  tidyr::complete(group, type, fill = list(count = 0)) %>%
  mutate(group = factor(group, levels = gene_type_levels),
         type  = factor(type,  levels = type_levels))

pD <- PlotNovelQTLCount(novel_qtl_count)

saveRDS(pD, file.path(OBJ_DIR, "D.rds"))
message("Panel D built.")
