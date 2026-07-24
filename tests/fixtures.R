# Builds a mock Fig3 project in a temp directory: a manifest plus tiny input
# files with the same schemas as the real ones. Sourcing scripts set the FIG3_*
# environment variables from the returned paths BEFORE sourcing fig3/00_setup.R.
#
# Three genes:
#   ENSG001 - fine-mapped in all three assays, colocalized in all three pairs
#   ENSG002 - fine-mapped in all three assays, colocalized in one pair
#   ENSG003 - eQTL only (so it is excluded from genes_all_three_assays)
#
# CREDIBLE-SET ID CONVENTION. A FastENLOC Signal is "<gene>:<cs_id>", and the
# pipeline recovers cs_id by stripping the "^ENSG[0-9]+:" prefix, then joins it
# back to the SuSiE cs_id column. So the SuSiE cs_id must be the BARE id, with
# no gene prefix, or the component join silently matches nothing. The pQTL cs_id
# additionally carries a "<UniProt>_" prefix, which is where the gene/protein map
# comes from.
make_fig3_fixture <- function(root = file.path(tempdir(), "fig3_fixture")) {
  unlink(root, recursive = TRUE)
  for (d in c("finemapping/coloc", "enloc", "annotated", "reference"))
    dir.create(file.path(root, d), recursive = TRUE, showWarnings = FALSE)

  wtsv <- function(df, p) write.table(df, p, sep = "\t", row.names = FALSE, quote = FALSE)

  # Nine variants, alternating between the two shared genes. Enough rows for the
  # per-MAF-bin regressions in panels C and D to actually fit.
  nv     <- 9
  pos    <- seq_len(nv) * 100
  snps   <- sprintf("chr1_%d_A_G", pos)
  genes  <- rep(c("ENSG001", "ENSG002"), length.out = nv)
  uni    <- c(ENSG001 = "P00001", ENSG002 = "P00002")

  # Bare credible-set ids, one per gene per assay.
  e_cs <- paste0("E_", genes, "_L1")
  s_cs <- paste0("S_", genes, "_L1")
  p_cs <- paste0(uni[genes], "_L1")

  # --- SuSiE fine-mapping ---------------------------------------------------
  eqtl <- data.frame(
    molecular_trait_id = paste0(genes, ".1"),        # version suffix is stripped
    variant        = snps,
    pip            = seq(0.9, 0.1, length.out = nv),
    cs_index       = "L1",
    posterior_mean = seq(0.8, 0.1, length.out = nv),
    posterior_sd   = 0.1,
    cs_id          = e_cs,
    gene_name      = paste0("GENE_", genes)
  )
  # An eQTL-only gene, so genes_all_three_assays is a strict subset.
  eqtl <- rbind(eqtl, data.frame(
    molecular_trait_id = "ENSG003.1", variant = "chr3_100_G_A", pip = 0.6,
    cs_index = "L1", posterior_mean = 0.3, posterior_sd = 0.1,
    cs_id = "E_ENSG003_L1", gene_name = "GENE_ENSG003"))
  wtsv(eqtl, file.path(root, "finemapping/coloc/eQTL_susie.COMB.annotated.cleaned.tsv.gz"))

  sqtl <- data.frame(
    molecular_trait_id = paste0(genes, ":clu_1"),    # gene is extracted from this
    variant        = snps,
    pip            = seq(0.8, 0.2, length.out = nv),
    posterior_mean = seq(-0.4, 0.4, length.out = nv),   # signed, to exercise abs()
    cs_index       = "L1",
    posterior_sd   = 0.1,
    cs_id          = s_cs
  )
  wtsv(sqtl, file.path(root, "finemapping/coloc/sQTL_susie.COMB.annotated.cleaned.tsv.gz"))

  pqtl <- data.frame(
    molecular_trait_id = genes,
    variant        = snps,
    pip            = seq(0.9, 0.3, length.out = nv),
    cs_index       = "L1",
    posterior_mean = seq(0.7, 0.2, length.out = nv),
    posterior_sd   = 0.1,
    cs_id          = p_cs,
    group          = "COMB"
  )
  # A non-COMB row that filter(group == "COMB") must drop.
  pqtl <- rbind(pqtl, transform(pqtl[1, ], group = "AFR"))
  wtsv(pqtl, file.path(root, "finemapping/coloc/pQTL_susie.COMB.annotated.cleaned.tsv.gz"))

  # --- FastENLOC gene level (a normal TSV) ----------------------------------
  gene_out <- function(p, glcp) wtsv(
    data.frame(Gene = c("ENSG001", "ENSG002"), GLCP = glcp), p)
  gene_out(file.path(root, "enloc/eQTL_pQTL.enloc.gene.out"), c(0.99, 0.98))
  gene_out(file.path(root, "enloc/sQTL_pQTL.enloc.gene.out"), c(0.99, 0.01))
  gene_out(file.path(root, "enloc/eQTL_sQTL.enloc.gene.out"), c(0.99, 0.01))

  # --- FastENLOC signal level (whitespace-aligned; the header is discarded) --
  # Signal = "<gene>:<left cs_id>(@)<gene>:<right cs_id>".
  sig_out <- function(p, left, right, rcp) writeLines(c(
    "Signal Num_SNP CPIP_qtl CPIP_gwas_marginal CPIP_gwas_qtl_prior RCP LCP",
    sprintf("%s(@)%s 5 0.9 0.1 0.8 %s 0.9", left, right, rcp)), p)
  g2 <- c("ENSG001", "ENSG002")
  sig_out(file.path(root, "enloc/eQTL_pQTL.enloc.sig.out"),
          paste0(g2, ":E_", g2, "_L1"), paste0(g2, ":", uni[g2], "_L1"),
          c("0.99", "0.98"))
  sig_out(file.path(root, "enloc/sQTL_pQTL.enloc.sig.out"),
          paste0(g2, ":S_", g2, "_L1"), paste0(g2, ":", uni[g2], "_L1"),
          c("0.99", "0.01"))
  sig_out(file.path(root, "enloc/eQTL_sQTL.enloc.sig.out"),
          paste0(g2, ":E_", g2, "_L1"), paste0(g2, ":S_", g2, "_L1"),
          c("0.99", "0.01"))

  # --- FastENLOC variant level ---------------------------------------------
  snp_out <- function(p, cs) writeLines(c(
    "Signal SNP PIP_qtl PIP_gwas_marginal PIP_gwas_qtl_prior SCP",
    sprintf("%s:%s %s 0.9 0.1 0.8 0.99", genes, cs, snps)), p)
  snp_out(file.path(root, "enloc/eQTL_pQTL.enloc.snp.out"), e_cs)
  snp_out(file.path(root, "enloc/sQTL_pQTL.enloc.snp.out"), s_cs)
  snp_out(file.path(root, "enloc/eQTL_sQTL.enloc.snp.out"), e_cs)

  # --- AoU annotation table -------------------------------------------------
  # Three variants per MAF bin so each panel-D facet can fit a slope.
  # Variant 1 has gvs_eur_af = 0.90, which must FOLD to 0.10, leaving max_af at
  # 0.30 from gvs_afr_af — the assertion in tests/test_fig3_builders.R.
  afr <- c(0.30, 0.20, 0.10, 0.030, 0.020, 0.015, 0.005, 0.003, 0.002)
  eur <- afr; eur[1] <- 0.90
  ann <- data.frame(
    chrom = "chr1", pos = pos, ref = "A", alt = "G",
    gvs_afr_af = afr,
    gvs_amr_af = afr * 0.9,
    gvs_eas_af = afr * 0.8,
    gvs_eur_af = eur,
    gvs_mid_af = afr * 0.7,
    gvs_sas_af = afr * 0.6,
    gvs_oth_af = afr * 0.5,
    gvs_all_af = afr,
    gvs_max_af = pmax(afr, eur)
  )
  wtsv(ann, file.path(root, "annotated/AoU.COMB.AC10.AN95.biallelic.v9.annotations.tsv.bgz"))

  # --- Manifest -------------------------------------------------------------
  manifest <- file.path(root, "inputs.tsv")
  rows <- read.delim("Fig3/inputs.tsv", colClasses = "character")
  rows <- rows[!rows$key %in% c("protein_lod", "median_tpm", "genebayes_shet"), ]
  write.table(rows, manifest, sep = "\t", row.names = FALSE, quote = FALSE)

  list(root = root, manifest = manifest)
}

# Point fig3/00_setup.R at the fixture. Call BEFORE sourcing it.
use_fig3_fixture <- function(fx) {
  Sys.setenv(FIG3_INPUT    = fx$root,
             FIG3_MANIFEST = fx$manifest,
             FIG3_DATA     = file.path(fx$root, "data"),
             FIG3_OBJ      = file.path(fx$root, "objects"),
             FIG3_OUTDIR   = fx$root)
}

# --- Fig4 mock project -------------------------------------------------------
# Small TSV stand-ins for Fig4's tabular inputs, enough to exercise the pure
# transforms. The ClinGen parquet and gene-set derivation are tested from an
# in-memory frame (see tests/test_fig4_helpers.R). Task 8 extends this fixture
# with the r2_twas / conditional_twas / drug / ClinGen inputs for the full run.
make_fig4_fixture <- function(root = file.path(tempdir(), "fig4_fixture")) {
  unlink(root, recursive = TRUE)
  for (d in c("high_level_coloc_outputs/harmonized",
              "high_level_coloc_outputs/final_summaries", "reference"))
    dir.create(file.path(root, d), recursive = TRUE, showWarnings = FALSE)
  wtsv <- function(df, p) write.table(df, p, sep = "\t", row.names = FALSE, quote = FALSE)

  # coloc gene table: gene / trait / trait_category / any_coloc
  genes <- expand.grid(g = 1:8, trait = c("eosinophil", "LDL"),
                       stringsAsFactors = FALSE)
  gene_tbl <- data.frame(
    gene = paste0("ENSG", sprintf("%08d", genes$g)),
    trait = genes$trait,
    trait_category = ifelse(genes$trait == "eosinophil", "hematological", "lipid"),
    any_coloc = TRUE)
  # fwrite gzips by the .gz extension, so read_any()'s `gzip -dc` path works.
  data.table::fwrite(
    gene_tbl,
    file.path(root, "high_level_coloc_outputs/harmonized/harmonized_coloc.gene.tsv.gz"),
    sep = "\t")

  # coloc rate: one union row per trait
  rate_tbl <- data.frame(
    trait = c("eosinophil", "LDL"), layer = "union",
    pct_stringent_union = c(30, 20), pct_lenient_union = c(45, 33))
  wtsv(rate_tbl, file.path(root, "high_level_coloc_outputs/final_summaries/coloc_rate_by_trait.tsv"))

  # R^2 TWAS: two AF thresholds per phenotype so the pivot has `0` and `0.01`
  r2 <- do.call(rbind, lapply(1:12, function(i) data.frame(
    phenotype_id = paste0("g", i), AF_threshold = c(0, 0.01),
    r.squared = c(0.30 + i/100, 0.20 + i/100), H2 = 0.5, SE = 0.05, Pval = 1e-8)))
  wtsv(r2, file.path(root, "R2_TWAS_combined.tsv"))

  # Conditional TWAS: enough colocalized + non-colocalized rows for MatchIt and glm
  # Conditional TWAS: 200 genes, seeded for reproducibility. any_coloc and the
  # gene-set memberships below are independent random draws at moderate prevalence
  # so that (a) Panel C's MatchIt has both a treated and a control group, and
  # (b) Panel D's per-outcome enrichment glms are well-posed (no constant outcome
  # or perfect separation, which would break broom's confint profiling).
  set.seed(2026)
  ng      <- 200
  gene_id <- paste0("ENSG", sprintf("%08d", 1:ng))
  desc    <- paste0("GENE", 1:ng)
  ct <- data.frame(
    status = "ok", gene = gene_id, coloc_trait = "eosinophil",
    molecular_trait_id = paste0("m", 1:ng),
    any_coloc = sample(c(TRUE, FALSE), ng, replace = TRUE),
    TWAS_path = "p1",
    z_rare_cond = stats::rnorm(ng, 2), rare_zscore = stats::rnorm(ng, 2),
    common_zscore = stats::rnorm(ng, 1), Vcor = 0.1, GLCP = 0.9,
    p_rare_cond = stats::runif(ng), rare_pvalue = stats::runif(ng),
    common_pvalue = stats::runif(ng), full_pvalue = stats::runif(ng),
    mean_TPM = stats::runif(ng, 1, 50), shet = stats::runif(ng, 0.001, 0.1),
    n_rare = sample(1:10, ng, replace = TRUE), Description = desc)
  # data.table::fwrite handles the .gz extension by gzipping
  data.table::fwrite(ct, file.path(root, "ConditionalTWAS.cleaned.tsv.gz"), sep = "\t")

  # Drug interactions: ~55% druggable, ~40% of those approved (independent)
  is_drug <- sample(c(TRUE, FALSE), ng, replace = TRUE, prob = c(.55, .45))
  drug_genes <- desc[is_drug]
  wtsv(data.frame(gene_name = drug_genes,
                  approved = sample(c(TRUE, FALSE), length(drug_genes),
                                    replace = TRUE, prob = c(.4, .6))),
       file.path(root, "reference/interactions.tsv"))

  # ClinGen parquet: ~65% Definitive/Strong (Mendelian), ~35% anemia/immune (heme)
  if (requireNamespace("arrow", quietly = TRUE)) {
    dir.create(file.path(root, "reference/evidence_clingen"), showWarnings = FALSE)
    clingen <- data.frame(
      targetFromSourceId = paste0("S", 1:ng),
      allelicRequirements = "autosomal recessive",
      confidence = sample(c("Definitive", "Strong", "Limited"), ng,
                          replace = TRUE, prob = c(.35, .30, .35)),
      score = 1, studyId = paste0("st", 1:ng), targetId = gene_id,
      diseaseFromSource = sample(c("hereditary anemia", "cardiomyopathy",
                                   "immune deficiency"), ng,
                                 replace = TRUE, prob = c(.35, .40, .25)))
    arrow::write_parquet(clingen,
      file.path(root, "reference/evidence_clingen/part-0.parquet"))
  }

  list(root = root)
}

use_fig4_fixture <- function(fx) {
  Sys.setenv(FIG4_INPUT    = fx$root,
             FIG4_DATA     = file.path(fx$root, "data"),
             FIG4_OBJ      = file.path(fx$root, "objects"),
             FIG4_OUTDIR   = fx$root)
}

# --- SuppFig1 mock project ---------------------------------------------------
# Mock AoU / GTEx GCT TPM matrices (Name, Description, then sample columns) and a
# small GENCODE gene-type table, enough to exercise the pure transforms. The GTF
# import and enrichR calls are not run here (rtracklayer/enrichR/network).
make_suppfig1_fixture <- function(root = file.path(tempdir(), "suppfig1_fixture")) {
  unlink(root, recursive = TRUE)
  for (d in c("expression", "external", "reference"))
    dir.create(file.path(root, d), recursive = TRUE, showWarnings = FALSE)
  wtsv <- function(df, p) write.table(df, p, sep = "\t", row.names = FALSE, quote = FALSE)

  # Six genes spanning the gene types the panels filter on, shared by both
  # cohorts, plus a 7th gene present only in AoU (absent from the GTEx GCT) so
  # join_mean_tpms's inner join has something cohort-exclusive to drop.
  ids  <- c("ENSG001.1","ENSG002.1","ENSG003.1","ENSG004.1","ENSG005.1","ENSG006.1")
  desc <- c("GENEA","TFGENE","LNC1","SNR1","HBB","RRNA1")
  types <- c("protein_coding","protein_coding","lncRNA","snRNA","protein_coding","rRNA")

  aou_only_id   <- "ENSG007.1"
  aou_only_desc <- "AOUONLY7"
  aou_only_type <- "protein_coding"

  aou <- data.frame(Name = c(ids, aou_only_id), Description = c(desc, aou_only_desc),
                    s1 = c(10, 8, 0.5, 0.3, 40, 0.05, 5),
                    s2 = c(12, 9, 0.4, 0.2, 44, 0.05, 7))
  wtsv(aou, file.path(root, "expression/aou.gct.tsv"))   # plain tsv stand-in

  gtex <- data.frame(Name = ids, Description = desc,
                     s1 = c(11, 2, 0.02, 0.25, 38, 0.06),
                     s2 = c(9,  3, 0.03, 0.30, 42, 0.05))
  wtsv(gtex, file.path(root, "external/gtex.gct.tsv"))

  ann <- data.frame(gene_id = c(ids, aou_only_id), gene_type = c(types, aou_only_type))
  saveRDS(ann, file.path(root, "reference/gene_annotations.rds"))

  list(root = root, aou_tpms = file.path(root, "expression/aou.gct.tsv"),
       gtex_tpms = file.path(root, "external/gtex.gct.tsv"),
       annotations = ann)
}
