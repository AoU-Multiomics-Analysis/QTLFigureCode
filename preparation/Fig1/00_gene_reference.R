# Extract stable gene IDs and annotations from the GENCODE v48 GTF.
source("preparation/common.R")
source("preparation/Fig1/functions.R")
config <- prep_init()
relative <- "output/qtl_circos_simplified_options/gene_reference.tsv"
message("Figure 1: Read GENCODE v48 genes from the GTF.")
genes <- fig1_read_gene_gtf(prep_input("input", "gencode.v48.annotation.gtf.gz"))
stopifnot(nrow(genes) == 78686L)
compared <- prep_has_reference("paper", relative)
if (compared) {
  saved <- read_tsv(prep_reference("paper", relative), show_col_types = FALSE)
  stopifnot(nrow(anti_join(genes, saved, by = names(genes))) == 0L,
            nrow(anti_join(saved, genes, by = names(genes))) == 0L)
}
prep_write_tsv(genes, "paper", relative)
message("Prepared 78,686 GENCODE gene records. Saved-reference comparison: ",
        if (compared) "passed" else "not run (reference absent)", ".")
