# Fig3 symmetric colocalization background + panel E rare-variant enrichment.
# Run from the project root:  Rscript tests/test_fig3_symmetric.R
#
# Everything here runs on synthetic frames: the symmetric universe needs the
# controlled-access SuSiE tables and a GENCODE GTF, neither of which is present
# locally, so the I/O builders are out of scope and the parsing / fitting /
# plotting seams are tested instead.
source("tests/fixtures.R")
source("tests/helpers.R")

fx <- make_fig3_fixture(file.path(tempdir(), "fig3_symmetric"))
use_fig3_fixture(fx)
source("fig3/00_setup.R")

# --- enloc Signal parsing ----------------------------------------------------
# Both traps below returned an EMPTY join rather than an error when wrong, which
# is why they are pinned.

ok("strips the gene prefix the SuSiE cs_id lacks (expression)",
   identical(.strip_gene_prefix("ENSG00000186827:ENSG00000186827.12_L4"),
             "ENSG00000186827.12_L4"))
ok("strips it from a UniProt-prefixed proteomics id",
   identical(.strip_gene_prefix("ENSG00000179403:Q6PCB0_ENSG00000179403.12_L2"),
             "Q6PCB0_ENSG00000179403.12_L2"))
ok("strips it from a splicing id carrying an intron coordinate block",
   identical(.strip_gene_prefix(
     "ENSG00000116288:chr1:7961793:7962763:clu_352_+:ENSG00000116288.13_L5"),
     "chr1:7961793:7962763:clu_352_+:ENSG00000116288.13_L5"))
ok("leaves an already-bare cs_id unchanged",
   identical(.strip_gene_prefix("ENSG00000000419.14_L2"), "ENSG00000000419.14_L2"))

# ONE pQTL credible set: the internal "_ENSG" belongs to the id itself.
# Splitting on a plain "_ENSG" returns c("Q6PCB0", "ENSG...L2") and matches
# nothing downstream.
ok("a single UniProt-prefixed partner id is NOT shredded at its internal _ENSG",
   identical(.partner_cs_ids("ENSG00000179403:Q6PCB0_ENSG00000179403.12_L2"),
             "Q6PCB0_ENSG00000179403.12_L2"))

two <- paste0("ENSG00000157870:Q8TBF2_ENSG00000157870.17_L5",
              "_ENSG00000157873:Q92956_ENSG00000157873.18_L2")
ok("two partner credible sets split on '_' before a NEW gene-prefixed id",
   identical(.partner_cs_ids(two),
             c("Q8TBF2_ENSG00000157870.17_L5", "Q92956_ENSG00000157873.18_L2")))

ok("a signal with no partner half yields no ids (the common, non-error case)",
   length(.partner_cs_ids(c("", ""))) == 0L)
ok("no signals at all yields no ids",
   length(.partner_cs_ids(character(0))) == 0L)

ok("gene accession survives the splicing table's strand prefix",
   identical(.bare_gene("-:ENSG00000000419"), "ENSG00000000419"))
ok("a bare accession passes through",
   identical(.bare_gene("ENSG00000111647"), "ENSG00000111647"))
ok("a field with no accession gives NA rather than a wrong gene",
   is.na(.bare_gene("no accession here")))

# --- fit_maf_enrichment ------------------------------------------------------
maf_frame <- function(n = 400, effect = 0, seed = 1) {
  set.seed(seed)
  x <- rep(c(TRUE, FALSE), each = n / 2)
  data.frame(
    marked  = x,
    null    = rep(c(TRUE, FALSE), times = n / 2),
    # effect > 0 shifts marked variants to LOWER log10 MAF, i.e. rarer
    log_maf = rnorm(n, mean = -1 - effect * x, sd = 0.3),
    layers  = factor(rep("eQTL", n)),
    pip     = rep(0.5, n))
}

r_eff <- fit_maf_enrichment(maf_frame(effect = 0.5), "marked")
ok("one row per fitted annotation", nrow(r_eff) == 1L)
ok("a rarer-variant effect comes back above 1", r_eff$fold > 1)
ok("its confidence interval excludes no-effect", r_eff$lo > 1)
ok("the magnitude is recovered", abs(r_eff$fold - 10^0.5) < 0.1)

r_null <- fit_maf_enrichment(maf_frame(effect = 0), "null")
ok("a null annotation spans 1", r_null$lo < 1 && r_null$hi > 1)

# The outcome is log10 MAF, so the fold change is 10^coef. exp() understates
# every effect -- a real bug in the exploratory code (see Fig3/NOTES.md).
d10 <- maf_frame(effect = 0)
d10$log_maf <- -1 - log10(2) * d10$marked          # exactly a 2-fold shift
r10 <- fit_maf_enrichment(d10, "marked")
ok("the fold change is base 10, matching the log10 outcome",
   abs(r10$fold - 2) < 1e-4)
ok("it is NOT the base-e back-transform",
   abs(r10$fold - exp(log10(2))) > 0.1)

dfloor <- maf_frame(effect = 0.5)
dfloor$rareflag <- c(rep(TRUE, 3), rep(FALSE, nrow(dfloor) - 3))
ok("an annotation below the variant floor is dropped",
   nrow(fit_maf_enrichment(dfloor, "rareflag", min_variants = 5)) == 0L)
ok("and fitted once the floor allows it",
   nrow(fit_maf_enrichment(dfloor, "rareflag", min_variants = 3)) == 1L)

# --- collapse_to_variants ----------------------------------------------------
# ~31% of colocalizing variants sit in BOTH layers' credible sets; without the
# collapse they enter the model twice as though independent.
u <- data.frame(
  variant = c("v1", "v1", "v2"),
  layer   = c("eQTL", "pQTL", "eQTL"),
  pip     = c(0.04, 0.64, 0.25),
  log_maf = c(-1, -1, -2),
  marked  = c(TRUE, FALSE, TRUE))
v <- collapse_to_variants(u, "marked")
ok("a two-layer variant is counted once", nrow(v) == 2L)
ok("its weight is the GEOMETRIC mean of the two PIPs",
   abs(v$pip[v$variant == "v1"] - sqrt(0.04 * 0.64)) < 1e-12)
ok("the layers it was fine-mapped in are recorded",
   v$layers[v$variant == "v1"] == "eQTL+pQTL")
ok("a single-layer variant records just its own layer",
   v$layers[v$variant == "v2"] == "eQTL")
ok("an annotation TRUE in either layer is TRUE for the variant",
   isTRUE(v$marked[v$variant == "v1"]))

# --- build_panelE_maf --------------------------------------------------------
enr <- function(...) data.frame(..., stringsAsFactors = FALSE)

e_sig <- enr(feature = c("dELS", "conserved"), fold = c(0.8, 2),
             lo = c(0.7, 1.8), hi = c(0.9, 2.3), n = c(100L, 100L),
             padj = c(1e-4, 0.9))
p_sig <- build_panelE_maf(e_sig)
ok("build_panelE_maf returns a ggplot", inherits(p_sig, "ggplot"))
ok("non-significant annotations are dropped", nrow(p_sig$data) == 1L)
ok("the surviving annotation is the significant one",
   as.character(p_sig$data$feature) == "dELS")

e_thin <- enr(feature = "conserved", fold = 2, lo = 1.8, hi = 2.3,
              n = 3L, padj = 1e-4)
ok("an annotation below the variant floor is not plotted",
   identical(build_panelE_maf(e_thin, min_n = 5), blank_plot))

# A CI end outside the scale is dropped by geom_linerange with only a warning,
# so the limits are derived from the data rather than hard-coded.
e_wide <- enr(feature = c("stop_gained", "dELS"), fold = c(6.6, 0.8),
              lo = c(4.5, 0.76), hi = c(9.7, 0.85), n = c(8L, 2000L),
              padj = c(1e-20, 1e-16))
# Assert the PROPERTY, not the internal representation: build the plot and
# confirm no linerange row was dropped for falling outside the scale.
b_wide <- ggplot2::ggplot_build(build_panelE_maf(e_wide))
lr <- b_wide$data[[3]]          # layer 3 = geom_linerange
ok("no confidence interval is dropped by the x limits",
   nrow(lr) == nrow(e_wide) && !any(is.na(lr$xmin)) && !any(is.na(lr$xmax)))

# An annotation missing from PANEL_E_LABELS would vanish silently -- exactly how
# stop_gained/frameshift disappeared during development.
e_unlab <- enr(feature = "not_a_known_annotation", fold = 2, lo = 1.8,
               hi = 2.3, n = 50L, padj = 1e-4)
msgs <- character(0)
withCallingHandlers(build_panelE_maf(e_unlab),
                    message = function(m) msgs <<- c(msgs, conditionMessage(m)))
ok("a significant but unlabelled annotation is named, not silently dropped",
   any(grepl("significant but unlabelled", msgs)))

e_empty <- enr(feature = character(0), fold = numeric(0), lo = numeric(0),
               hi = numeric(0), n = integer(0), padj = numeric(0))
ok("an empty enrichment table renders the blank placeholder",
   identical(build_panelE_maf(e_empty), blank_plot))

cat("test_fig3_symmetric.R: all assertions passed\n")
