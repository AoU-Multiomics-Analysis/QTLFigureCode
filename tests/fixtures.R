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
