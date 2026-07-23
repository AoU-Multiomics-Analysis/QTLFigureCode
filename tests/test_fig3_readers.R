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

# Declared type policy for the CS-level reader: Signal is the only identifier
# column; everything else (Num_SNP, the CPIP_* columns, RCP, LCP) is numeric.
ok("normal path: Signal is character",    is.character(d$Signal))
ok("normal path: Num_SNP is numeric",     is.numeric(d$Num_SNP))
ok("normal path: CPIP_qtl is numeric",    is.numeric(d$CPIP_qtl))
ok("normal path: CPIP_gwas_marginal is numeric", is.numeric(d$CPIP_gwas_marginal))
ok("normal path: CPIP_gwas_qtl_prior is numeric", is.numeric(d$CPIP_gwas_qtl_prior))
ok("normal path: RCP is numeric",         is.numeric(d$RCP))
ok("normal path: LCP is numeric",         is.numeric(d$LCP))

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

# --- Finding 1: exact threshold round-trip, boundary < 0.5 -------------------
# Four rows with posteriors 0.9, 0.3, 0.3, 0.05. Sorted descending, local FDRs
# are 0.10, 0.70, 0.70, 0.95; cumulative means are 0.10, 0.40, 0.50, 0.6125.
# At a 50% control level, the first three rows qualify (0.50 <= 0.5), so the
# boundary is the tied pair of 0.3 rows: cutoff = 0.3. In IEEE-754,
# `1 - (1 - 0.3)` is 0.30000000000000004, NOT bitwise-equal to 0.3, so the old
# `thresh <- 1 - lfdr[k]` reconstruction makes `probs >= thresh` wrongly
# exclude BOTH 0.3 rows (0.3 >= 0.30000000000000004 is FALSE). The fix reads
# the cutoff back from `probs` itself, so the boundary rows are included.
boundary_df <- tibble::tibble(p = c(0.9, 0.3, 0.3, 0.05))
boundary_lab <- add_fdr_threshold_label(boundary_df, "p", control_level = 0.5, quiet = TRUE)
ok("high posterior above the boundary is labelled TRUE", isTRUE(boundary_lab$Coloc[1]))
ok("boundary probability below 0.5 is labelled TRUE (row 2)", isTRUE(boundary_lab$Coloc[2]))
ok("boundary probability below 0.5 is labelled TRUE (row 3)", isTRUE(boundary_lab$Coloc[3]))
ok("tie at the cutoff: both tied rows labelled the same",
   boundary_lab$Coloc[2] == boundary_lab$Coloc[3])
ok("row below the cutoff is labelled FALSE (row 4)", isFALSE(boundary_lab$Coloc[4]))

# --- Finding 1: NA probability is labelled FALSE, not NA ---------------------
na_df <- tibble::tibble(p = c(0.99, NA, 0.01))
na_lab <- add_fdr_threshold_label(na_df, "p", control_level = 0.05, quiet = TRUE)
ok("NA probability is labelled FALSE (not NA)", isFALSE(na_lab$Coloc[2]))
ok("Coloc column has no NAs", !anyNA(na_lab$Coloc))

# --- Finding 2: a ragged row (wrong token count) errors, not silently corrupts
ragged <- file.path(tmp, "pair.enloc.ragged.out")
writeLines(c(
  "Signal  Num_SNP  CPIP_qtl  CPIP_gwas_marginal  CPIP_gwas_qtl_prior  RCP  LCP",
  "ENSG001:chr1_100(@)ENSG001:chr1_200   5  0.90  0.10  0.80  0.95  0.90",
  "ENSG002:chr2_100(@)ENSG002:chr2_200   3  0.50  0.20  0.40  0.10"   # missing a column
), ragged)
expect_error("a ragged row errors instead of silently recycling",
             LoadEnlocTable(ragged, ENLOC_CS_COLS))

# --- Finding 3: a header-only file returns 0 rows with the right column names
header_only <- file.path(tmp, "pair.enloc.empty.out")
writeLines(c(
  "Signal  Num_SNP  CPIP_qtl  CPIP_gwas_marginal  CPIP_gwas_qtl_prior  RCP  LCP"
), header_only)
empty <- LoadEnlocTable(header_only, ENLOC_CS_COLS)
ok("header-only file returns 0 rows",           nrow(empty) == 0L)
ok("header-only file has the right col names",  identical(names(empty), ENLOC_CS_COLS))

# The blind spot that let the original bug through: names/row count alone
# don't catch a header-only result coming back all-character. Assert the
# actual column types match the same identifier/measurement policy as the
# normal (non-empty) path above.
ok("header-only file: Signal is character",      is.character(empty$Signal))
ok("header-only file: Num_SNP is numeric",       is.numeric(empty$Num_SNP))
ok("header-only file: CPIP_qtl is numeric",      is.numeric(empty$CPIP_qtl))
ok("header-only file: CPIP_gwas_marginal is numeric", is.numeric(empty$CPIP_gwas_marginal))
ok("header-only file: CPIP_gwas_qtl_prior is numeric", is.numeric(empty$CPIP_gwas_qtl_prior))
ok("header-only file: RCP is numeric",           is.numeric(empty$RCP))
ok("header-only file: LCP is numeric",           is.numeric(empty$LCP))

# The actual consequence described in the finding: load_input() row-binds
# every file for a manifest key with purrr::map_dfr(), so a header-only
# FastENLOC file and a populated one for the same key must bind_rows()
# cleanly, with the empty file contributing zero rows and no type conflict.
combined <- dplyr::bind_rows(empty, d)
ok("bind_rows(header-only, normal) succeeds and has the normal row count",
   nrow(combined) == nrow(d))
ok("bind_rows result: Signal is character",  is.character(combined$Signal))
ok("bind_rows result: Num_SNP is numeric",    is.numeric(combined$Num_SNP))
ok("bind_rows result: RCP is numeric",        is.numeric(combined$RCP))
ok("bind_rows result: LCP is numeric",        is.numeric(combined$LCP))

# Also verify the variant-level (.snp.out) column policy on the header-only
# path: Signal and SNP are identifiers (character), the rest are numeric.
snp_header_only <- file.path(tmp, "pair.enloc.snp.empty.out")
writeLines(c(
  "Signal  SNP  PIP_qtl  PIP_gwas_marginal  PIP_gwas_qtl_prior  SCP"
), snp_header_only)
empty_v <- LoadEnlocTable(snp_header_only, ENLOC_VARIANT_COLS)
ok("variant header-only: Signal is character", is.character(empty_v$Signal))
ok("variant header-only: SNP is character",    is.character(empty_v$SNP))
ok("variant header-only: PIP_qtl is numeric",  is.numeric(empty_v$PIP_qtl))
ok("variant header-only: SCP is numeric",      is.numeric(empty_v$SCP))
combined_v <- dplyr::bind_rows(empty_v, v)
ok("bind_rows(variant header-only, normal) succeeds",
   nrow(combined_v) == nrow(v))

cat("test_fig3_readers.R: all assertions passed\n")
