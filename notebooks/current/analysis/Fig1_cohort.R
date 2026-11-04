# %% Count post-QC samples and retain pre-QC aggregate availability
#| The participant-availability export is absent. Pre-QC availability starts from its supplied aggregate table.
#| RNA totals are counted from the sample list. Protein totals are counted from the matrix and covariate headers.
folder <- "output/figure1_flowchart_draft/science_print_size"
pre_qc <- read_tsv(prep_input("paper", file.path(folder, "cohort_counts_2026-10-26.tsv")),
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
cohort_totals <- tibble(category = c("Total", "Post-QC"),
                RNA = c(before$n[before$name == "has_rna"], nrow(samples)),
                Proteome = c(before$n[before$name == "has_proteome"], length(protein_ids)))

cohort <- pre_qc %>% transmute(assay = recode(name, has_rna = "RNA", has_proteome = "Proteome"),
  ancestry = recode(str_to_upper(ancestry_pred_other), OTH = "REM"), n) %>%
  mutate(ancestry = factor(ancestry, levels = ancestry_levels))
message("Counted 8,895 RNA and 9,581 protein samples after QC. No participant records are exported.")
