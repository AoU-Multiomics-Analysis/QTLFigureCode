# Fig1 alternative enrichment — glm(annotated ~ PIP) over the SuSiE tables.
# Run from the project root:  Rscript tests/test_fig1_enrichment_glm.R
#
# The fixture reproduces the shape of the real tables, including the trap that
# broke the first run: in the splicing table tss / gene_name / distTSS /
# Constrained / lof.pLI are NA for every row, so fread types them logical and a
# `where(is.logical)` column selection collects them as annotations -- one of
# which is all-NA after the pivot and leaves glm() nothing to fit.
source("tests/helpers.R")
suppressMessages(source("fig1/00_setup.R"))
source("fig1/enrichment_glm_fit.R")

set.seed(1)
n <- 600

# PIP-dependent annotation: high-PIP variants carry `missense_variant` far more
# often, so the slope must come out clearly positive. `intron_variant` runs the
# other way, and `synonymous_variant` is independent of PIP.
pip <- rep(seq(0.02, 0.98, length.out = n / 2), 2)
missense   <- runif(n) < pip * 0.8
intronic   <- runif(n) < (1 - pip) * 0.8
synonymous <- runif(n) < 0.3

consequence <- vapply(seq_len(n), function(i) {
  terms <- c(if (missense[i]) "missense_variant",
             if (intronic[i]) "intron_variant",
             if (synonymous[i]) "synonymous_variant")
  if (!length(terms)) NA_character_ else paste(terms, collapse = ",")
}, character(1))

fixture <- data.frame(
  variant     = paste0("chr1_", seq_len(n), "_A_G"),
  pip         = pip,
  consequence = consequence,
  # The seven real flag columns.
  CTCF.only     = runif(n) < 0.1,
  CTCF.bound    = runif(n) < 0.1,
  pELS          = runif(n) < 0.2,
  DNase.H3K4me3 = runif(n) < 0.1,
  dELS          = runif(n) < 0.2,
  PLS           = runif(n) < pip * 0.5,
  FANTOM5       = runif(n) < 0.15,
  # The all-NA columns that the splicing table really has. fread types every one
  # of these as logical.
  tss         = NA,
  gene_name   = NA,
  distTSS     = NA,
  Constrained = NA,
  lof.pLI     = NA,
  stringsAsFactors = FALSE
)

path <- file.path(tempdir(), "glm_fixture.tsv")
data.table::fwrite(fixture, path, sep = "\t")

# The fixture must actually reproduce the trap, or the test below proves nothing.
typed <- data.table::fread(path)
ok("the fixture reproduces the all-NA-typed-logical trap",
   all(vapply(typed[, c("tss", "gene_name", "distTSS", "Constrained", "lof.pLI")],
              is.logical, logical(1))))

fits <- suppressWarnings(GlmAnnotationEnrichment(path))
slopes <- fits[fits$term == "pip", ]

ok("the run completes despite the all-NA columns", nrow(slopes) > 0)
ok("no all-NA column is fitted as an annotation",
   !any(c("tss", "gene_name", "distTSS", "Constrained", "lof.pLI") %in% slopes$name))
ok("every cCRE flag is fitted", all(CCRE_FLAG_COLS %in% slopes$name))
ok("VEP consequences are one-hot expanded out of the comma-separated string",
   all(c("missense_variant", "intron_variant", "synonymous_variant") %in% slopes$name))

get_slope <- function(nm) slopes$estimate[slopes$name == nm]
ok("a PIP-driven annotation gets a positive slope", get_slope("missense_variant") > 0)
ok("a PIP-depleted annotation gets a negative slope", get_slope("intron_variant") < 0)
ok("a PIP-independent annotation is not significant",
   slopes$p.value[slopes$name == "synonymous_variant"] > 0.05)
ok("results are sorted by descending estimate",
   !is.unsorted(rev(fits$estimate)))

# A renamed or missing flag must be an error, not an annotation that vanishes.
dropped <- file.path(tempdir(), "glm_fixture_dropped.tsv")
data.table::fwrite(fixture[, setdiff(names(fixture), "PLS")], dropped, sep = "\t")
expect_error("a missing flag column is an error, not a silent drop",
             GlmAnnotationEnrichment(dropped))

retyped <- fixture
retyped$FANTOM5 <- ifelse(retyped$FANTOM5, "yes", "no")
rt_path <- file.path(tempdir(), "glm_fixture_retyped.tsv")
data.table::fwrite(retyped, rt_path, sep = "\t")
expect_error("a flag column that stops being logical is an error",
             GlmAnnotationEnrichment(rt_path))

cat("test_fig1_enrichment_glm.R: all assertions passed\n")
