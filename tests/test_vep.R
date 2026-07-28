# Behaviour of the shared VEP consequence expansion (common/vep.R).
# Run from the project root:  Rscript tests/test_vep.R
suppressPackageStartupMessages({ library(dplyr) })
source("common/vep.R")
source("tests/helpers.R")

d <- data.frame(
  variant     = paste0("v", 1:4),
  consequence = c("intron_variant,NMD_transcript_variant",
                  "missense_variant",
                  ".",
                  "splice_region_variant,intron_variant"),
  PLS         = c(TRUE, FALSE, TRUE, FALSE),
  stringsAsFactors = FALSE)
e <- expand_vep_consequences(d)

ok("adds one column per distinct term",
   all(c("intron_variant", "NMD_transcript_variant", "missense_variant",
         "splice_region_variant") %in% names(e)))
ok("added columns are logical",
   all(vapply(e[c("intron_variant", "missense_variant")], is.logical, logical(1))))
ok("a multi-term row sets EVERY one of its terms",
   isTRUE(e$splice_region_variant[4]) && isTRUE(e$intron_variant[4]))
ok("a term absent from a row is FALSE, not NA",
   isFALSE(e$missense_variant[1]))
ok("VEP's '.' placeholder produces no column", !"." %in% names(e))
ok("an existing column is never clobbered", identical(e$PLS, d$PLS))

# A frame with no consequence column passes straight through.
ok("no consequence column is a no-op",
   identical(expand_vep_consequences(data.frame(x = 1)), data.frame(x = 1)))
ok("zero rows is a no-op",
   nrow(expand_vep_consequences(d[0, ])) == 0L)

cat("test_vep.R: all assertions passed\n")
