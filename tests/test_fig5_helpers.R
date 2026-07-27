# Fig5 pure-seam unit tests. Run from the project root (via tests/run_all.R).
# Sources fig5/00_setup.R for the folded seams; uses in-memory fixture frames so
# no controlled-access data, Excel, or gzip reads are needed.
source("fig5/00_setup.R")
source("tests/fixtures.R")
source("tests/helpers.R")

fx <- make_fig5_fixture()

# --- derive_clingen ---------------------------------------------------------
cg <- derive_clingen(fx$clingen_raw, fx$susie_gene_map)
ok("derive_clingen keeps only Definitive/Strong AD/AR/SD genes",
   setequal(unique(cg$targetFromSourceId), c("HEXA", "MLC1", "GENE_AD", "GENE_BOTH")))
ok("derive_clingen drops Limited-confidence and non-AD/AR/SD rows",
   !any(c("LIMGENE", "XLGENE") %in% cg$targetFromSourceId))
ok("derive_clingen keeps a gene curated under BOTH AD and AR as two rows",
   sum(cg$targetFromSourceId == "GENE_BOTH") == 2L &&
     setequal(cg$allelicRequirements[cg$targetFromSourceId == "GENE_BOTH"], c("AD", "AR")))
ok("derive_clingen maps gene_name -> Ensembl targetId via the SuSiE map",
   cg$targetId[cg$targetFromSourceId == "HEXA"][1] == "ENSGH")

# --- annotate_eqtl ----------------------------------------------------------
clingen_genes <- cg$targetFromSourceId
eqtl <- annotate_eqtl(fx$susie, fx$genebayes, clingen_genes)
ok("annotate_eqtl flags ClinGen genes",
   all(eqtl$ClinGen[eqtl$gene_name %in% c("HEXA", "MLC1", "GENE_AD", "GENE_BOTH")]))
ok("annotate_eqtl leaves non-ClinGen genes unflagged",
   !any(eqtl$ClinGen[eqtl$gene_name %in% c("CTRL1", "NONCLIN")]))
ok("annotate_eqtl sets rare = gvs_max_af < 0.01",
   identical(eqtl$rare, eqtl$gvs_max_af < 0.01) && !eqtl$rare[eqtl$gvs_max_af == 0.02][1])
ok("annotate_eqtl joins GeneBayes post_mean (CTRL1 is low-shet)",
   abs(eqtl$post_mean[eqtl$gene_name == "CTRL1"][1] - 0.005) < 1e-9)

# --- build_variant_freq -----------------------------------------------------
vf <- build_variant_freq(eqtl)
ok("build_variant_freq returns one row per distinct variant", nrow(vf) == 7L)
ok("build_variant_freq normalizes the variant key (colons -> underscores)",
   all(!grepl(":", vf$variant_key)))
ok("build_variant_freq prefixes rare variants with 'rare'",
   all(grepl("^rare", vf$freq_label[vf$rare])))
ok("build_variant_freq formats the aFC label with a sign",
   any(grepl("aFC -6", vf$afc_label)))

# --- id utilities -----------------------------------------------------------
ok("make_variant_id normalizes delimiters and upcases alleles",
   make_variant_id("chr1", 100, "a", "g") == "chr1_100_A_G")
ok("normalize_variant_id collapses repeated separators",
   normalize_variant_id("chr1::100**A") == "chr1_100_A")

cat("test_fig5_helpers.R passed.\n")
