suppressPackageStartupMessages(library(tidyverse))
stopifnot(file.exists("notebooks/current/analysis/S14_helpers.R"))
source("notebooks/current/analysis/S14_helpers.R")
notebook_paths <- c("notebooks/current/supplements/S14.R", "notebooks/current/analysis/S14_helpers.R", "notebooks/current/plots/S14_discovery.R")
notebook_code <- paste(unlist(lapply(notebook_paths, readLines, warn = FALSE)), collapse = "\n")
stopifnot(!grepl("readRDS|saveRDS|prep_output|prep_save_rds|parse\\(|eval\\(|source\\(", notebook_code))

gtf <- tempfile(fileext = ".gtf")
writeLines(c("# GENCODE test", 'chr1\ttest\tgene\t10\t20\t.\t+\t.\tgene_id "ENSG000001.2"; gene_name "A"; gene_type "protein_coding";',
  'chr1\ttest\texon\t10\t15\t.\t+\t.\tgene_id "ENSG000001.2";',
  'chr2\ttest\tgene\t30\t40\t.\t-\t.\tgene_id "ENSG000002.1"; gene_name "B"; gene_type "lncRNA";'), gtf)
genes <- s14_read_gene_gtf(gtf)
stopifnot(identical(genes$gene_id, c("ENSG000001", "ENSG000002")))
raw <- tibble(molecular_trait_id = c("ENSG000001.2", "ENSG000001.2", "ENSG000001.2", "ENSG000002.1"),
  cs_id = c("A_L1", "A_L1", "A_L2", "B_L1"))
aou <- discovery_sets(raw, "molecular_trait_id", "AoU", "COMB", genes)
interval <- discovery_sets(raw[1, ], "molecular_trait_id", "INTERVAL", "study", genes)
topmed <- discovery_sets(raw[1, ], "molecular_trait_id", "TOPMed", "study", genes)
stopifnot(nrow(aou) == 3L, sum(aou$protein_coding) == 2L,
  discovery_counts(aou)$credible_sets == 3L,
  discovery_counts(aou)$genes_with_credible_sets == 2L,
  identical(shared_discovery_genes(bind_rows(aou, interval, topmed)), "ENSG000001"))
cloud_error <- tryCatch(s14_read_gene_gtf("gs://bucket/input.gtf"), error = conditionMessage)
stopifnot(grepl("localization", cloud_error))
tissues <- tibble(tissue = LETTERS[1:6], sample_size = c(100, 200, 150, 250, 180, 130),
  credible_sets = c(1000, 1300, 1600, 2000, 1400, 1100))
envelope <- discovery_envelope(tissues)
stopifnot(nrow(envelope$boundary) == 721L, all(envelope$membership$normalized_distance2 < 1),
  !envelope$details$confidence_region)
unlink(gtf)
message("Inline S14 source-count, GTF and ellipse checks passed.")
