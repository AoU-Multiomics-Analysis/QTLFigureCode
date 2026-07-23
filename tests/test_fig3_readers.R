# Tests for Fig3's FastENLOC readers and FDR labelling.
# Run from the project root:  Rscript tests/test_fig3_readers.R
source("fig3/00_setup.R")
source("tests/helpers.R")

tmp <- file.path(tempdir(), "fig3_readers")
unlink(tmp, recursive = TRUE); dir.create(tmp, recursive = TRUE)

# enloc .sig.out: a header line the parser discards, then whitespace-aligned rows.
sig <- file.path(tmp, "pair.enloc.sig.out")
writeLines(c(
  "Signal  Num_SNP  CPIP_qtl  CPIP_gwas_marginal  CPIP_gwas_qtl_prior  RCP  LCP",
  "ENSG001:chr1_100(@)ENSG001:chr1_200   5  0.90  0.10  0.80  0.95  0.90",
  "ENSG002:chr2_100(@)ENSG002:chr2_200   3  0.50  0.20  0.40  0.10  0.05"
), sig)

d <- LoadEnlocTable(sig, ENLOC_CS_COLS)
ok("reader drops the header line",        nrow(d) == 2L)
ok("reader applies the column names",     identical(names(d), ENLOC_CS_COLS))
ok("numeric columns are converted",       is.numeric(d$RCP))
ok("integer columns are converted",       is.numeric(d$Num_SNP))
ok("Signal stays character",              is.character(d$Signal))
ok("values are parsed positionally",      d$RCP[1] == 0.95)

# A .snp.out has six columns, not seven.
snp <- file.path(tmp, "pair.enloc.snp.out")
writeLines(c(
  "Signal  SNP  PIP_qtl  PIP_gwas_marginal  PIP_gwas_qtl_prior  SCP",
  "ENSG001:chr1_100  chr1_100_A_G  0.9  0.1  0.8  0.99",
  "ENSG002:chr2_100  chr2_100_C_T  0.2  0.1  0.1  0.01"
), snp)

v <- LoadEnlocTable(snp, ENLOC_VARIANT_COLS)
ok("variant-level reader parses 6 columns", identical(names(v), ENLOC_VARIANT_COLS))
ok("variant-level SCP is numeric",          is.numeric(v$SCP))

expect_error("wrong column count errors", LoadEnlocTable(snp, ENLOC_CS_COLS))

# FDR labelling: with one near-certain and one near-zero posterior, only the
# former should be labelled at a 5% Bayesian FDR.
lab <- add_fdr_threshold_label(d, "RCP", quiet = TRUE)
ok("adds the Coloc column",              "Coloc" %in% names(lab))
ok("Coloc is logical",                   is.logical(lab$Coloc))
ok("high posterior is labelled TRUE",    isTRUE(lab$Coloc[1]))
ok("low posterior is labelled FALSE",    isFALSE(lab$Coloc[2]))

lab2 <- add_fdr_threshold_label(d, "RCP", label_col = "Sig", quiet = TRUE)
ok("label_col is honoured",              "Sig" %in% names(lab2))

expect_error("an unknown column errors", add_fdr_threshold_label(d, "NOPE", quiet = TRUE))

cat("test_fig3_readers.R: all assertions passed\n")
