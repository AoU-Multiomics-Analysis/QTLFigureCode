#!/usr/bin/env Rscript
# =============================================================================
# Regenerate Fig5's gene-model reference tables from a GENCODE GTF.
#
#   Fig5/reference/gene_spans.tsv   one row per gene  (the gene record)
#   Fig5/reference/gene_exons.tsv   exons, merged across transcripts
#
# Both are headerless 5-column TSVs read by Fig5/00_setup.R as
#   gene, chrom, start, end, strand
# and drawn by gene_model_track(): spans become the gene line + TSS arrow
# (so a gene MUST appear exactly once), exons become the filled blocks.
#
# Exons are reduced per gene because GENCODE lists them once per transcript;
# without merging, shared exons are drawn dozens of times on top of each other.
#
# Coordinates stay 1-based inclusive, exactly as GENCODE reports them.
#
# Usage (from the project root):
#   Rscript tools/make_gene_models.R --gtf path/to/gencode.v48.annotation.gtf.gz
#   Rscript tools/make_gene_models.R --gtf <gtf> --genes ACVRL1,ANKRD33,SPR
#   Rscript tools/make_gene_models.R --gtf <gtf> --check   # verify, write nothing
#
# These tables are GENCODE-derived and contain no participant data, which is why
# they are committed to the repo.
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })

GENES  <- c("ACVRL1", "ANKRD33", "SPR")   # the genes Fig5's panels D-G display
OUTDIR <- "Fig5/reference"

# --- args -------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop("missing value after ", flag)
  args[i + 1L]
}
gtf     <- arg("--gtf", Sys.getenv("GENCODE_GTF", ""))
genes   <- strsplit(arg("--genes", paste(GENES, collapse = ",")), ",")[[1]]
outdir  <- arg("--out", OUTDIR)
check   <- "--check" %in% args

if (!nzchar(gtf))
  stop("no GTF given. Pass --gtf <path> or set GENCODE_GTF.\n",
       "  GENCODE: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/")
if (!file.exists(gtf)) stop("GTF not found: ", gtf)

message("GTF   : ", gtf)
message("genes : ", paste(genes, collapse = ", "))

# --- extract ----------------------------------------------------------------
source("common/gene_models.R")
gm    <- extract_gene_models(gtf, genes = genes, gene_types = NULL)
spans <- gm$spans
exons <- gm$exons

message("spans : ", nrow(spans), " rows")
message("exons : ", nrow(exons), " rows (merged across transcripts)")

# --- write, or check against what is committed ------------------------------
targets <- list(gene_spans = spans, gene_exons = exons)

for (nm in names(targets)) {
  path <- file.path(outdir, paste0(nm, ".tsv"))
  if (check) {
    if (!file.exists(path)) { message("CHECK ", nm, ": no committed file at ", path); next }
    old <- fread(path, col.names = c("gene", "chrom", "start", "end", "strand"),
                 colClasses = "character")
    new <- targets[[nm]][, lapply(.SD, as.character)]
    same <- isTRUE(all.equal(as.data.frame(old), as.data.frame(new)))
    message("CHECK ", nm, ": ", if (same) "identical to the committed file"
            else "DIFFERS from the committed file")
    if (!same) {
      message("  committed: ", nrow(old), " rows | regenerated: ", nrow(new), " rows")
      message("  (a different GENCODE release is the usual cause)")
    }
  } else {
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    fwrite(targets[[nm]], path, sep = "\t", col.names = FALSE)
    message("wrote ", path)
  }
}
