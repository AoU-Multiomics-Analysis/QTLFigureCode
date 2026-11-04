# Prepare the Figure 1A count inputs. Do not write participant identifiers.
# Pre-QC availability counts are a supplied aggregate source table.
# Post-QC totals are counted from the RNA list and protein matrix/covariate headers.
source("preparation/common.R")
config <- prep_init()
folder <- "output/figure1_flowchart_draft/science_print_size"
pre_qc <- read_tsv(prep_reference("paper", file.path(folder, "cohort_counts_2026-10-26.tsv")),
                   show_col_types = FALSE)
samples <- read_tsv(prep_input("input", "reference/COMB_eQTL_9k_sample_list.tsv"),
                   col_types = cols_only(research_id = col_character()), progress = FALSE)
bed_file <- prep_input_file("QTL_FIG1_PROTEIN_BED_FILE", "input", "reference/AoU_comb_10750_pQTL.bed.gz")
covariate_file <- prep_input_file("QTL_FIG1_PROTEIN_COVARIATES_FILE", "input",
                                "reference/AoU_COMB_10750_pQTL_covars_no_genetic_pcs.tsv")
bed_connection <- gzfile(bed_file, "rt")
bed_header <- str_split(readLines(bed_connection, n = 1L), "\t")[[1]]
close(bed_connection)
covariate_header <- str_split(readLines(covariate_file, n = 1L), "\t")[[1]]
stopifnot(identical(bed_header[1:4], c("#chr", "start", "end", "gene_id")), covariate_header[[1]] == "ID")
protein_ids <- bed_header[-(1:4)]
stopifnot(!anyDuplicated(samples$research_id), !anyDuplicated(protein_ids),
          setequal(protein_ids, covariate_header[-1]), nrow(samples) == 8895L, length(protein_ids) == 9581L)
before <- pre_qc %>% group_by(name) %>% summarise(n = sum(n), .groups = "drop")
stopifnot(before$n[before$name == "has_rna"] == 9130L,
          before$n[before$name == "has_proteome"] == 10056L)
totals <- tibble(category = c("Total", "Post-QC"),
                RNA = c(before$n[before$name == "has_rna"], nrow(samples)),
                Proteome = c(before$n[before$name == "has_proteome"], length(protein_ids)))
saved <- read_tsv(prep_reference("paper", file.path(folder,
  "Figure1_Science_184mm_seven_categories_postQC_20261027_cohort_totals.tsv")), show_col_types = FALSE)
stopifnot(isTRUE(all.equal(totals, saved, check.attributes = FALSE)))
prep_write_tsv(pre_qc, "paper", file.path(folder, "cohort_counts_2026-10-26.tsv"))
prep_write_tsv(totals, "paper", file.path(folder,
  "Figure1_Science_184mm_seven_categories_postQC_20261027_cohort_totals.tsv"))
message("Figure 1A verified: RNA 8,895 and proteome 9,581 post-QC samples; identifiers were not written.")
