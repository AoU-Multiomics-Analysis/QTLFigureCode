suppressPackageStartupMessages({library(tidyverse); library(patchwork)})
source("notebooks/current/plots/Fig5_locus.R")
models <- list(spans = tibble(gene = c("A", "B", "C"), chrom = c("chr12", "chr2", "chr12"), start = c(90, 90, 300), end = c(120, 120, 350), strand = c("+", "+", "-")),
  exons = tibble(gene = c("A", "A", "B", "C"), chrom = c("chr12", "chr12", "chr2", "chr12"), start = c(90, 115, 90, 300), end = c(95, 120, 95, 305), strand = "+"))
selected <- gene_models_for_window("chr12", 100, 150, models)
stopifnot(identical(selected$spans$gene, "A"), identical(selected$exons$start, 115))
stopifnot(identical(regional_window_bp(c(100, 200, 300), c(FALSE, TRUE, FALSE), 200, 20, 5), c(180, 220)))
stopifnot(nrow(gene_models_for_window("chr12", 400, 500, models)$spans) == 0L)
stopifnot(identical(names(cs_palette(c("cs1", "cs2", "cs1"))), c("cs1", "cs2")))
gene_spans <- models$spans
gene_exons <- models$exons
qtl <- tibble(position = c(100, 200, 300), pip = c(0.2, 0.95, 0.3), cs_id = "cs1", rare = TRUE)
locus <- fig5_plot_acvrl1_locus(qtl, tibble(start = numeric(), end = numeric(), type = character()), focal_position = 200)
stopifnot(inherits(patchwork::patchworkGrob(locus), "gtable"))
gtf <- tempfile(fileext = ".gtf")
writeLines(c(
  'chr12\tGENCODE\tgene\t90\t150\t.\t+\t.\tgene_name "A"; gene_type "protein_coding";',
  'chr12\tGENCODE\texon\t90\t110\t.\t+\t.\tgene_name "A"; gene_type "protein_coding";',
  'chr12\tGENCODE\texon\t100\t120\t.\t+\t.\tgene_name "A"; gene_type "protein_coding";',
  'chr12\tGENCODE\texon\t121\t130\t.\t+\t.\tgene_name "A"; gene_type "protein_coding";',
  'chr2\tGENCODE\tgene\t90\t150\t.\t+\t.\tgene_name "B"; gene_type "protein_coding";',
  'chr12\tGENCODE\tgene\t90\t150\t.\t+\t.\tgene_name "C"; gene_type "lncRNA";'), gtf)
fresh_models <- s31_read_gene_models(gtf, "chr12", 100, 140)
stopifnot(identical(fresh_models$spans$gene, "A"),
          nrow(fresh_models$exons) == 1L, fresh_models$exons$start == 90L,
          fresh_models$exons$end == 130L)
unlink(gtf)
message("Inline Figure 5 locus checks passed.")
