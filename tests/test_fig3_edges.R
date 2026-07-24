# Regression tests for the two defects in Fig3.R that broke the credible-set graph
# behind panel B.
#
#   Fig3.R:462  str_extract(sig_1, "^EN[<43;65;42MSG[0-9]+") - a pasted terminal
#               escape sequence inside the regex. Under stringr/ICU this pattern
#               does not compile at all: ICU reads the inner [0-9] as a nested
#               character class, leaving the outer [ unclosed, so str_extract
#               raises U_REGEX_MISSING_CLOSE_BRACKET and the script halts there.
#               ([<43;65;42M is an SGR mouse-click report, so it was most likely
#               pasted in by a stray terminal event after the analysis last ran.)
#   Fig3.R:477  bind_rows() listed cleaned_eqtl_pqtl twice and never bound
#               cleaned_eqtl_sqtl, so eQTL-sQTL edges were silently absent. Of the
#               two, this is the one that could have reached published output.
#
# coloc_cs_edges() is a pure function, so no manifest or fixture is needed.
# Run from the project root:  Rscript tests/test_fig3_edges.R
source("fig3/00_setup.R")
source("tests/helpers.R")

cs_colocs <- tibble::tribble(
  ~type,       ~Signal,                                                ~cs_1,      ~cs_2,      ~RCP, ~LCP, ~Coloc,
  "eQTL-pQTL", "ENSG001:chr1_100(@)ENSG001:chr1_200",                  "chr1_100", "chr1_200", 0.90, 0.90, TRUE,
  "eQTL-sQTL", "ENSG002:chr2_100(@)ENSG002:chr2_200",                  "chr2_100", "chr2_200", 0.80, 0.80, TRUE,
  "sQTL-pQTL", "ENSG003:chr3_100(@)ENSG003:chr3_200",                  "chr3_100", "chr3_200", 0.70, 0.70, TRUE,
  "eQTL-pQTL", "ENSG004:chr4_100(@)ENSG005:chr4_200",                  "chr4_100", "chr4_200", 0.60, 0.60, TRUE,
  "eQTL-pQTL", "ENSG006:chr6_100(@)ENSG006:chr6_200",                  "chr6_100", "chr6_200", 0.50, 0.50, FALSE,
  "eQTL-sQTL", "ENSG007:chr7_100(@)ENSG007:chr7_200_ENSG007:chr7_300", "chr7_100", NA,         0.95, 0.95, TRUE
)

e <- coloc_cs_edges(cs_colocs)

# Defect 2: every assay pair must contribute, exactly once.
ok("all three assay pairs are represented",
   setequal(e$type, c("eQTL-pQTL", "eQTL-sQTL", "sQTL-pQTL")))
ok("eQTL-sQTL is present (was dropped by the duplicated bind_rows)",
   "ENSG002" %in% e$gene)
ok("no pair is bound twice", sum(e$gene == "ENSG003") == 1L)

# Defect 1: the eQTL-pQTL gene regex must match.
ok("eQTL-pQTL edges survive the gene regex", "ENSG001" %in% e$gene)

# Behaviour that must be preserved.
ok("signals spanning two genes are dropped",  !"ENSG004" %in% e$gene)
ok("non-colocalized signals are dropped",     !"ENSG006" %in% e$gene)
ok("a multi-signal second side is unnested",  sum(e$gene == "ENSG007") == 2L)
ok("unnested credible sets lose the gene prefix",
   setequal(e$cs_2[e$gene == "ENSG007"], c("chr7_200", "chr7_300")))
ok("edges carry gene/type/cs columns",
   identical(sort(names(e)), sort(c("gene", "RCP", "LCP", "type", "cs_1", "cs_2"))))
ok("total edge count is as expected", nrow(e) == 5L)

cat("test_fig3_edges.R: all assertions passed\n")
