# Resolve a manifest key to the path that manifest DECLARES for it, creating the
# parent directory. Fixtures must follow the manifest rather than hardcode a
# layout: the two silently drifted apart the first time the real manifests were
# repointed at a different staging root, and every affected test failed with a
# "Missing local input" that pointed at the manifest, not at the fixture.
manifest_path <- function(rows, root, key, type = NULL) {
  r <- rows[rows$key == key, , drop = FALSE]
  if (!is.null(type)) r <- r[r$type == type, , drop = FALSE]
  if (nrow(r) != 1L)
    stop("fixture: key '", key, "'", if (!is.null(type)) paste0(" / type '", type, "'") else "",
         " matched ", nrow(r), " manifest rows; expected exactly 1")
  out <- file.path(root, r$local_subpath)
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  out
}

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
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  rows <- read.delim("Fig3/inputs.tsv", colClasses = "character")
  rows <- rows[!rows$key %in% c("protein_lod", "median_tpm", "genebayes_shet"), ]
  mp <- function(key, type = NULL) manifest_path(rows, root, key, type)

  wtsv <- function(df, p) write.table(df, p, sep = "\t", row.names = FALSE, quote = FALSE)

  # --- minimal GENCODE GTF --------------------------------------------------
  # fig3/symmetric.R needs Ensembl-canonical CDS records to decide whether a
  # coding call refers to the eGene's own dominant protein. One canonical CDS
  # per gene, covering the fixture's variant positions (chr1:100-900).
  gtf_path <- mp("gene_gtf")
  dir.create(dirname(gtf_path), recursive = TRUE, showWarnings = FALSE)
  gtf_line <- function(gene, start, end, canonical = TRUE)
    sprintf(paste0("chr1\tHAVANA\tCDS\t%d\t%d\t.\t+\t0\t",
                   "gene_id \"%s.1\"; transcript_id \"%s_T1\";%s\n"),
            start, end, gene, gene,
            if (canonical) " tag \"Ensembl_canonical\";" else "")
  gz <- gzfile(gtf_path, "w")
  writeLines(c("##description: fixture", "##provider: none"), gz)
  cat(gtf_line("ENSG001", 50, 550),  file = gz)
  cat(gtf_line("ENSG002", 560, 950), file = gz)
  # a NON-canonical CDS overlapping ENSG001's span: must never satisfy the
  # canonical-CDS test on its own
  cat(gtf_line("ENSG003", 100, 200, canonical = FALSE), file = gz)
  close(gz)

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

  # --- Per-ancestry allele frequencies --------------------------------------
  # These ride along on the annotated SuSiE tables (Fig3 no longer reads a
  # separate annotation table). Three variants per MAF bin so each panel-D
  # facet can fit a slope. Variant 1 has gvs_eur_af = 0.90, which must FOLD to
  # 0.10, leaving max_af at 0.30 from gvs_afr_af — the assertion in
  # tests/test_fig3_builders.R.
  afr <- c(0.30, 0.20, 0.10, 0.030, 0.020, 0.015, 0.005, 0.003, 0.002)
  eur <- afr; eur[1] <- 0.90
  af_cols <- data.frame(
    gvs_afr_af = afr,
    gvs_amr_af = afr * 0.9,
    gvs_eas_af = afr * 0.8,
    gvs_eur_af = eur,
    gvs_mid_af = afr * 0.7,
    gvs_sas_af = afr * 0.6,
    gvs_oth_af = afr * 0.5
  )
  # Every assay reports the SAME frequencies for a given variant; a disagreement
  # is what .variant_af_map()'s conflict guard exists to catch.
  af_for <- function(v) af_cols[match(v, snps), , drop = FALSE]

  # Annotation columns for the enrichment feature set. Deterministic from the
  # variant index so assertions can target specific rows. Variant 1 is PLS +
  # missense + conserved + TSS-proximal; variant 2 carries two VEP terms.
  annot_for <- function(v) {
    i <- seq_along(v)
    data.frame(
      PLS             = i %% 2 == 1,
      pELS            = i %% 3 == 0,
      dELS            = i %% 4 == 0,
      `CTCF.bound`    = i %% 5 == 0,
      `CTCF.only`     = FALSE,
      `DNase.H3K4me3` = i %% 6 == 0,
      FANTOM5         = i %% 7 == 0,
      consequence     = ifelse(i == 1, "missense_variant",
                        ifelse(i == 2, "intron_variant,NMD_transcript_variant",
                               "intron_variant")),
      phylop          = ifelse(i %% 2 == 1, 3.5, 0.2),
      distTSS         = ifelse(i %% 2 == 1, 500L, 50000L),
      gene_type       = ifelse(i %% 3 == 0, "lncRNA", "protein_coding"),
      cs_size         = 2L + (i %% 4),
      # Read by fig3/symmetric.R. gene_id is the ENSG the credible set belongs
      # to (strand-prefixed in the real splicing table, hence .bare_gene()).
      # SpliceAI scores gate the "clean canonical" coding redefinition: only
      # variant 1 is splice-capable, so it is the one dropped from
      # missense_variant by the SPLICEAI_MAX criterion.
      splice_ai_acceptor_gain_score = ifelse(i == 1, 0.8, 0.01),
      splice_ai_acceptor_loss_score = 0.01,
      splice_ai_donor_gain_score    = 0.01,
      splice_ai_donor_loss_score    = 0.01,
      check.names     = FALSE, stringsAsFactors = FALSE)
  }

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
  eqtl <- cbind(eqtl, af_for(eqtl$variant))
  # annot_for() is index-based (not variant-matched like af_for()), so it must
  # be applied here -- before the ENSG003 row below -- and NOT recomputed
  # after further rows are appended, or a variant repeated at a different row
  # position (e.g. sqtl_extra below) would pick up a different annotation and
  # spuriously trip the conflicting-values guard.
  eqtl <- cbind(eqtl, annot_for(eqtl$variant))
  # An eQTL-only gene, so genes_all_three_assays is a strict subset. Its variant
  # appears in no other table and carries NA frequencies, exercising the
  # all-NA drop in .variant_af_map().
  eqtl <- rbind(eqtl, data.frame(
    molecular_trait_id = "ENSG003.1", variant = "chr3_100_G_A", pip = 0.6,
    cs_index = "L1", posterior_mean = 0.3, posterior_sd = 0.1,
    cs_id = "E_ENSG003_L1", gene_name = "GENE_ENSG003",
    af_cols[NA_integer_, , drop = FALSE],
    annot_for(1L)[1, ]))

  # A variant fine-mapped for TWO different genes with DIFFERENT gene-relative
  # annotations -- the exact real-data shape (33,435 variants over 100,000
  # sampled) that broke .variant_feature_map()'s original all-columns guard.
  # Reuses variant 5 (snps[5] = "chr1_500_A_G", already a credible-set member
  # of ENSG001 above with distTSS = 500 / gene_type = "protein_coding, i.e.
  # TSS-proximal + protein-coding -- deliberately NOT variant 1/snps[1], which
  # a separate test below (test_fig3_builders.R's "same variant, different
  # genes" case) independently duplicates under ENSG002 with distTSS/gene_type
  # left UNCHANGED; reusing snps[1] here too would make that duplication and
  # this one disagree about the same (variant, gene) pair and trip the very
  # guard both are meant to demonstrate does NOT fire spuriously), now ALSO
  # fine-mapped for ENSG002 with a DIFFERENT gene-relative annotation
  # (distTSS = 50000 / gene_type = "lncRNA", i.e. TSS-distal + not
  # protein-coding). Its variant-level values (AF, consequence, phylop) are
  # left IDENTICAL to what it already carries for ENSG001 -- af_for() looks
  # them up by variant, and dual_gene_annot only overrides distTSS/gene_type --
  # so the variant-level guard in .variant_feature_map() stays silent and only
  # the gene-relative map (.variant_gene_feature_map()) sees the two differing
  # rows.
  dual_gene_variant <- snps[5]
  # NOTE: annot_for() indexes by seq_along(v), so annot_for(5L) would give row
  # 1's annotation again (seq_along of a length-1 vector is always 1), not row
  # 5's -- pull row 5 straight out of the already-annotated `eqtl` frame
  # (still un-shifted at this point: the ENSG003 row above was appended, not
  # inserted, so row 5 is still variant 5's original annotation) instead.
  dual_gene_annot <- eqtl[5, names(annot_for(1L)), drop = FALSE]
  dual_gene_annot$distTSS   <- 50000L
  dual_gene_annot$gene_type <- "lncRNA"
  eqtl <- rbind(eqtl, data.frame(
    molecular_trait_id = "ENSG002.1", variant = dual_gene_variant, pip = 0.5,
    cs_index = "L2", posterior_mean = 0.4, posterior_sd = 0.1,
    cs_id = "E_ENSG002_L2", gene_name = "GENE_ENSG002",
    af_for(dual_gene_variant),
    dual_gene_annot))
  # gene_id rides on the real SuSiE exports and is read by fig3/symmetric.R.
  # Added at write time, after every rbind(), so it cannot desynchronise them.
  # The splicing table's copy is strand-prefixed, as in the real export.
  add_gene_id <- function(df, strand_prefix = FALSE) {
    g <- stringr::str_extract(df$molecular_trait_id, "ENSG[0-9]+")
    df$gene_id <- if (strand_prefix) paste0("-:", g) else g
    df
  }

  wtsv(add_gene_id(eqtl), mp("susie_eqtl"))

  # NOTE: splicing ids are INTRON-level. `sqtl_extra` below adds a SECOND intron
  # cluster for ENSG001 that shares a variant with the first, reproducing the
  # duplicate (variant, gene) key that appears once these ids collapse to a gene.
  sqtl <- data.frame(
    molecular_trait_id = paste0(genes, ":clu_1"),    # gene is extracted from this
    variant        = snps,
    pip            = seq(0.8, 0.2, length.out = nv),
    posterior_mean = seq(-0.4, 0.4, length.out = nv),   # signed, to exercise abs()
    cs_index       = "L1",
    posterior_sd   = 0.1,
    cs_id          = s_cs
  )
  sqtl <- cbind(sqtl, af_for(sqtl$variant))
  # As with eqtl above, annot_for() must be applied before sqtl_extra is
  # derived: sqtl_extra copies row 1 wholesale (including its already-computed
  # annotation values), so the shared variant keeps ONE consistent annotation
  # instead of picking up a second, index-shifted one after the rbind below.
  sqtl <- cbind(sqtl, annot_for(sqtl$variant))
  # A second intron cluster of ENSG001, sharing variant 1 with clu_1 but with a
  # LARGER |posterior_mean|, so the "keep the strongest effect" collapse has a
  # determinate winner to assert on.
  sqtl_extra <- sqtl[1, ]
  sqtl_extra$molecular_trait_id <- "ENSG001:clu_2"
  sqtl_extra$cs_id              <- "S_ENSG001_L2"
  sqtl_extra$cs_index           <- "L2"
  sqtl_extra$posterior_mean     <- 0.95            # strongest for this (variant, gene)
  sqtl <- rbind(sqtl, sqtl_extra)
  wtsv(add_gene_id(sqtl, strand_prefix = TRUE), mp("susie_sqtl"))

  pqtl <- data.frame(
    molecular_trait_id = genes,
    variant        = snps,
    pip            = seq(0.9, 0.3, length.out = nv),
    cs_index       = "L1",
    posterior_mean = seq(0.7, 0.2, length.out = nv),
    posterior_sd   = 0.1,
    cs_id          = p_cs,
    group          = "COMB.proteomics"   # real tables carry the per-assay value
  )
  pqtl <- cbind(pqtl, af_for(pqtl$variant))
  pqtl <- cbind(pqtl, annot_for(pqtl$variant))
  wtsv(add_gene_id(pqtl), mp("susie_pqtl"))

  # --- FastENLOC gene level (a normal TSV) ----------------------------------
  gene_out <- function(p, glcp) wtsv(
    data.frame(Gene = c("ENSG001", "ENSG002"), GLCP = glcp), p)
  gene_out(mp("enloc_gene", "eQTL-pQTL"), c(0.99, 0.98))
  gene_out(mp("enloc_gene", "sQTL-pQTL"), c(0.99, 0.01))
  gene_out(mp("enloc_gene", "eQTL-sQTL"), c(0.99, 0.01))

  # --- FastENLOC signal level (whitespace-aligned; the header is discarded) --
  # Signal = "<gene>:<left cs_id>(@)<gene>:<right cs_id>".
  sig_out <- function(p, left, right, rcp) writeLines(c(
    "Signal Num_SNP CPIP_qtl CPIP_gwas_marginal CPIP_gwas_qtl_prior RCP LCP",
    sprintf("%s(@)%s 5 0.9 0.1 0.8 %s 0.9", left, right, rcp)), p)
  g2 <- c("ENSG001", "ENSG002")
  sig_out(mp("enloc_cs", "eQTL-pQTL"),
          paste0(g2, ":E_", g2, "_L1"), paste0(g2, ":", uni[g2], "_L1"),
          c("0.99", "0.98"))
  sig_out(mp("enloc_cs", "sQTL-pQTL"),
          paste0(g2, ":S_", g2, "_L1"), paste0(g2, ":", uni[g2], "_L1"),
          c("0.99", "0.01"))
  sig_out(mp("enloc_cs", "eQTL-sQTL"),
          paste0(g2, ":E_", g2, "_L1"), paste0(g2, ":S_", g2, "_L1"),
          c("0.99", "0.01"))

  # --- FastENLOC variant level ---------------------------------------------
  # `extra` appends additional, hand-written rows (same whitespace-aligned
  # shape) after the per-variant ones -- used below to give the shared
  # dual_gene_variant a SECOND (SNP, Gene) row under ENSG002, so
  # coloc_variant_level()'s gene-relative join has something real to match
  # against end-to-end, not just .variant_gene_feature_map() in isolation.
  snp_out <- function(p, cs, extra = character()) writeLines(c(
    "Signal SNP PIP_qtl PIP_gwas_marginal PIP_gwas_qtl_prior SCP",
    sprintf("%s:%s %s 0.9 0.1 0.8 0.99", genes, cs, snps),
    extra), p)
  dual_gene_row <- sprintf("ENSG002:E_ENSG002_L2 %s 0.9 0.1 0.8 0.99", dual_gene_variant)
  snp_out(mp("enloc_variant", "eQTL-pQTL"), e_cs, extra = dual_gene_row)
  snp_out(mp("enloc_variant", "sQTL-pQTL"), s_cs)
  snp_out(mp("enloc_variant", "eQTL-sQTL"), e_cs)

  # --- Manifest -------------------------------------------------------------
  manifest <- file.path(root, "inputs.tsv")
  write.table(rows, manifest, sep = "\t", row.names = FALSE, quote = FALSE)

  list(root = root, manifest = manifest)
}

# Path of a manifest key inside an already-built fixture, for tests that poke at
# an input file directly. Same rule as manifest_path(): never hardcode a layout.
fixture_path <- function(fx, key, type = NULL)
  manifest_path(read.delim(fx$manifest, colClasses = "character"), fx$root, key, type)

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
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  # Fig4 has no fixture manifest of its own -- it runs against the real
  # Fig4/inputs.tsv -- so every file must land where THAT declares it.
  rows <- read.delim("Fig4/inputs.tsv", colClasses = "character")
  mp <- function(key) manifest_path(rows, root, key)
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
    mp("coloc_gene"),
    sep = "\t")

  # coloc rate: one union row per trait
  rate_tbl <- data.frame(
    trait = c("eosinophil", "LDL"), layer = "union",
    pct_stringent_union = c(30, 20), pct_lenient_union = c(45, 33))
  wtsv(rate_tbl, mp("coloc_rate"))

  # R^2 TWAS: two AF thresholds per phenotype so the pivot has `0` and `0.01`
  r2 <- do.call(rbind, lapply(1:12, function(i) data.frame(
    phenotype_id = paste0("g", i), AF_threshold = c(0, 0.01),
    r.squared = c(0.30 + i/100, 0.20 + i/100), H2 = 0.5, SE = 0.05, Pval = 1e-8)))
  wtsv(r2, mp("r2_twas"))

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
  data.table::fwrite(ct, mp("conditional_twas"), sep = "\t")

  # Drug interactions: ~55% druggable, ~40% of those approved (independent)
  is_drug <- sample(c(TRUE, FALSE), ng, replace = TRUE, prob = c(.55, .45))
  drug_genes <- desc[is_drug]
  wtsv(data.frame(gene_name = drug_genes,
                  approved = sample(c(TRUE, FALSE), length(drug_genes),
                                    replace = TRUE, prob = c(.4, .6))),
       mp("drug_interactions"))

  # ClinGen parquet: ~65% Definitive/Strong (Mendelian), ~35% anemia/immune (heme)
  if (requireNamespace("arrow", quietly = TRUE)) {
    clingen_dir <- mp("clingen")
    dir.create(clingen_dir, recursive = TRUE, showWarnings = FALSE)
    clingen <- data.frame(
      targetFromSourceId = paste0("S", 1:ng),
      allelicRequirements = "autosomal recessive",
      confidence = sample(c("Definitive", "Strong", "Limited"), ng,
                          replace = TRUE, prob = c(.35, .30, .35)),
      score = 1, studyId = paste0("st", 1:ng), targetId = gene_id,
      diseaseFromSource = sample(c("hereditary anemia", "cardiomyopathy",
                                   "immune deficiency"), ng,
                                 replace = TRUE, prob = c(.35, .40, .25)))
    arrow::write_parquet(clingen, file.path(clingen_dir, "part-0.parquet"))
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

# --- Fig5 mock project -------------------------------------------------------
# In-memory frames (and small on-disk stand-ins) matching the schemas the folded
# Fig5 loaders read, enough to exercise the pure seams: derive_clingen,
# annotate_eqtl, build_variant_freq, variant_effect_by_inheritance, pct_at_element.
#
# Six genes: HEXA / MLC1 (ClinGen AR), GENE_AD (AD), GENE_BOTH (curated under
# BOTH AD and AR -- pins the preserved fan-out double-count), CTRL1 (low-shet
# control, not ClinGen), NONCLIN (not ClinGen, high shet -> excluded from B/C).
make_fig5_fixture <- function(root = file.path(tempdir(), "fig5_fixture")) {
  unlink(root, recursive = TRUE)
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  wtsv <- function(df, p) write.table(df, p, sep = "\t", row.names = FALSE, quote = FALSE)

  # SuSiE eQTL table (one rare variant per gene, plus a common HEXA variant to
  # exercise the rarity flag). `ref` is present on purpose (the data-masking trap).
  susie <- data.frame(
    gene_id      = c("ENSGH.1","ENSGM.2","ENSGAD.1","ENSGB.1","ENSGC.1","ENSGN.1","ENSGH.1"),
    gene_name    = c("HEXA","MLC1","GENE_AD","GENE_BOTH","CTRL1","NONCLIN","HEXA"),
    variant      = c("chr15:72346:C:T","chr22:50100:G:A","chr9:130000:A:C",
                     "chr1:200000:T:G","chr2:30000:C:G","chr3:40000:G:T","chr15:72999:A:G"),
    gvs_max_af   = c(0.002, 0.005, 0.003, 0.004, 0.006, 0.008, 0.02),
    gvs_max_subpop = c("gvs_eur_af","gvs_afr_af","gvs_sas_af","gvs_eas_af",
                       "gvs_amr_af","gvs_eur_af","gvs_eur_af"),
    log2_aFC       = c(-6, -4, -3, -2, 0.5, 1.0, -1),
    posterior_mean = c(0.9, 0.7, 0.6, 0.5, 0.2, 0.3, 0.2),
    pip            = c(0.98, 0.95, 0.92, 0.95, 0.95, 0.95, 0.99),
    ref            = c("C","G","A","T","C","G","A"),
    clinvar_classification = c("pathogenic","uncertain significance",NA,"benign",NA,NA,NA),
    stringsAsFactors = FALSE)
  wtsv(susie, file.path(root, "susie.tsv"))

  susie_gene_map <- susie %>% distinct(gene_id, gene_name) %>%
    mutate(gene_id = sub("\\..*", "", gene_id))

  # GeneBayes shet: CTRL1 low (< 0.01 -> a control gene), the rest high.
  genebayes <- data.frame(
    ensg = c("ENSGH","ENSGM","ENSGAD","ENSGB","ENSGC","ENSGN"),
    post_mean = c(0.8, 0.7, 0.6, 0.5, 0.005, 0.5), stringsAsFactors = FALSE)
  saveRDS(genebayes, file.path(root, "genebayes.rds"))

  # ClinGen raw (as fread(skip=6, header=FALSE, select=1:7) yields it): 6 filler
  # lines, then gene / hgnc / disease / mondo / moi / sop / confidence rows.
  # GENE_BOTH is curated AD *and* AR; LIMGENE (Limited) and XLGENE (XL) are dropped.
  clingen_lines <- c(
    "CLINGEN GENE DISEASE VALIDITY,,,,,,", "banner1,,,,,,", "banner2,,,,,,",
    "banner3,,,,,,", "banner4,,,,,,", "header,,,,,,",
    "HEXA,HGNC1,d1,MONDO1,AR,SOP1,Definitive",
    "MLC1,HGNC2,d2,MONDO2,AR,SOP2,Strong",
    "GENE_AD,HGNC3,d3,MONDO3,AD,SOP3,Definitive",
    "GENE_BOTH,HGNC4,d4,MONDO4,AD,SOP4,Definitive",
    "GENE_BOTH,HGNC5,d5,MONDO5,AR,SOP5,Definitive",
    "LIMGENE,HGNC6,d6,MONDO6,AD,SOP6,Limited",
    "XLGENE,HGNC7,d7,MONDO7,XL,SOP7,Definitive")
  writeLines(clingen_lines, file.path(root, "clingen.csv"))
  clingen_raw <- data.table::fread(file.path(root, "clingen.csv"), skip = 6,
    header = FALSE, select = 1:7,
    col.names = c("gene","hgnc","disease","mondo","moi","sop","confidence"))

  # SPR browser tracks for pct_at_element: one mark, alt = half of ref -> -50%.
  pos <- seq(-300, 300, by = 8)
  refv <- 10 * exp(-(pos / 60)^2)
  spr_tracks <- bind_rows(
    tibble(mark = "CTCF", type = factor("reference", levels = c("reference","alternate")),
           relative_pos = pos, score = refv),
    tibble(mark = "CTCF", type = factor("alternate", levels = c("reference","alternate")),
           relative_pos = pos, score = refv * 0.5))

  list(root = root, susie = susie, susie_gene_map = susie_gene_map,
       genebayes = genebayes, clingen_raw = clingen_raw, spr_tracks = spr_tracks,
       susie_path = file.path(root, "susie.tsv"),
       genebayes_path = file.path(root, "genebayes.rds"),
       clingen_path = file.path(root, "clingen.csv"))
}
