suppressPackageStartupMessages(library(tidyverse))
source("preparation/S14/functions.R")
genes <- tibble(gene_id = c("ENSG1", "ENSG2", "ENSG3"),
                gene_name = c("A", "B", "C"),
                gene_type = c("protein_coding", "protein_coding", "lncRNA"))
raw <- tibble(molecular_trait_id = c("ENSG1.2", "ENSG1.2", "ENSG1.2", "ENSG3.1"),
              cs_id = c("L1", "L1", "L2", "L1"))
x <- discovery_sets(raw, "molecular_trait_id", "AoU", "COMB", genes)
stopifnot(nrow(x) == 3L, sum(x$protein_coding) == 2L)
counts <- discovery_counts(x)
stopifnot(counts$credible_sets == 3L, counts$genes_with_credible_sets == 2L)
three <- bind_rows(x,
  discovery_sets(tibble(molecular_trait_id = c("ENSG1", "ENSG2"), cs_id = c(1, 1)),
                 "molecular_trait_id", "INTERVAL", "QTD001000", genes),
  discovery_sets(tibble(molecular_trait_id = "ENSG1", cs_id = 1),
                 "molecular_trait_id", "TOPMed", "Whole_blood.maf0001.cs.txt", genes))
stopifnot(identical(shared_discovery_genes(three), "ENSG1"))
# An ancestry-specific analysis must not narrow the three-study intersection.
three <- bind_rows(three, discovery_sets(
  tibble(molecular_trait_id = "ENSG2", cs_id = 1), "molecular_trait_id", "AoU", "EAS", genes))
stopifnot(identical(shared_discovery_genes(three), "ENSG1"))
message("S14 preparation tests passed.")
