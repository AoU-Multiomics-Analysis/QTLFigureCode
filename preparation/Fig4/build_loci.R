source("preparation/common.R")
source("preparation/Fig4/methods.R")
source("preparation/Fig4/current_methods.R")

message("Prepare unique gene–trait QTL membership for S24")
membership_input <- prep_input_file("QTL_FIG4_MEMBERSHIP", "input", c(
  "high_level_coloc_outputs/final_summaries/colocalizing_genes_long.tsv",
  "Fig4Data/high_level_coloc_outputs/final_summaries/colocalizing_genes_long.tsv"))
upset <- fig4_upset_counts(read_any(membership_input))
stopifnot(upset$multilayer_n == 1617L)
prep_write_tsv(upset$intersections, "source", "Fig4/coloc_upset_intersection_counts.tsv")
prep_write_tsv(upset$membership %>% select(-qtl_combination), "source", "Fig4/data/coloc_upset_gene_trait_membership.tsv")

message("Join CD40 fine-mapping variants to participating colocalization signals")
signal_path <- prep_input_file("QTL_FIG4_SIGNAL", "input", c(
  "high_level_coloc_outputs/harmonized/harmonized_coloc.signal.tsv.gz",
  "Fig4Data/high_level_coloc_outputs/harmonized/harmonized_coloc.signal.tsv.gz"))
signal <- read_any(signal_path)
qtl_files <- c(eQTL = "COMB.expression_SusieMerged.annotated.ancestry_skew.tsv.gz",
               pQTL = "COMB.proteomics_SusieMerged.annotated.ancestry_skew.tsv.gz",
               sQTL = "COMB.splicing_SusieMerged..annotated.ancestry_skew.tsv.gz")
gene_id <- "ENSG00000101017"
variants <- imap_dfr(qtl_files, function(relative, assay) {
  alternatives <- unique(c(relative, sub("SusieMerged[.][.]annotated", "SusieMerged.annotated", relative)))
  path <- prep_input_any("susie", alternatives)
  data <- as_tibble(data.table::fread(path, select = c("position", "variant", "cs_id", "pip", "rsid", "gene_id"), showProgress = FALSE))
  data <- if (assay == "sQTL") filter(data, str_detect(.data$gene_id, fixed(.env$gene_id))) else filter(data, .data$gene_id == .env$gene_id)
  data %>% mutate(layer = assay, qtl_cs = cs_id)
})
cd40 <- fig4_cd40_data(signal, variants, gene_id)
stopifnot(nrow(cd40) > 0)
prep_write_tsv(cd40, "source", "Fig4/cd40_coloc_locus_plot_data.tsv")
message("S24/S25 derived tables prepared without rerunning colocalization")
