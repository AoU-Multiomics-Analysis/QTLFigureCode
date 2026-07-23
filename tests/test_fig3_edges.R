# Regression tests for the two data-corrupting defects in Fig3.R.
#
#   Fig3.R:462  str_extract(sig_1, "^EN[<43;65;42MSG[0-9]+") - a pasted terminal
#               escape sequence inside the regex. gene_1 was NA for every
#               eQTL-pQTL row, so filter(gene_1 == gene_2) dropped all of them.
#   Fig3.R:477  bind_rows() listed cleaned_eqtl_pqtl twice and never bound
#               cleaned_eqtl_sqtl.
#
# Together those meant the credible-set graph behind panel B was built from
# sQTL-pQTL edges only.
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
