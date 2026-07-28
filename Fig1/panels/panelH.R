# Panel H — fine-mapping functional enrichment by assay
# Inputs: annot (COMB annotated SuSiE tables for eQTL / sQTL / pQTL)
source("fig1/00_setup.R")

# One enrichment computation per annotated table, labelled by its assay/type.
annot <- input_files("annot")
QTLFinemappingEnrichments <- purrr::pmap_dfr(
  list(annot$path, annot$Assay),
  function(p, ty) fread(p) %>% ComputeFinemappingAnnotationEnrichment() %>% mutate(type = ty)
)

# Names must match VEP's own spelling exactly (3_prime_UTR_variant, not
# _utr_): these are matched by an equality join, so a casing difference
# silently drops the annotation rather than erroring.
variant_functional_groups <- tibble::tribble(
  ~name, ~functional_group,
  # UTR
  "3_prime_UTR_variant", "Regulatory Annotation",
  "5_prime_UTR_variant", "Regulatory Annotation",
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
  "mature_miRNA_variant", "Noncoding RNA",
  # Gene proximity / background
  "downstream_gene_variant", "Intronic/Intergenic",
  "upstream_gene_variant", "Intronic/Intergenic",
  "intron_variant", "Intronic/Intergenic",
  "intergenic_variant", "Intronic/Intergenic"
)

pH <- PlotFinemappingEnrichment(QTLFinemappingEnrichments, variant_functional_groups)

saveRDS(pH, file.path(OBJ_DIR, "H.rds"))
message("Panel H built.")
