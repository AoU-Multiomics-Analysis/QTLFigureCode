# Tests for common/manifest.R. Run from the project root:
#   Rscript tests/test_manifest.R
suppressPackageStartupMessages({
  library(data.table); library(tibble); library(purrr)
})
source("tests/helpers.R")

tmp <- file.path(tempdir(), "manifest_test")
unlink(tmp, recursive = TRUE)
dir.create(file.path(tmp, "sub"), recursive = TRUE)

wtsv <- function(df, p) write.table(df, p, sep = "\t", row.names = FALSE, quote = FALSE)
wtsv(data.frame(v = 1:2), file.path(tmp, "sub/a.tsv"))
wtsv(data.frame(v = 3:5), file.path(tmp, "sub/b.tsv"))
wtsv(data.frame(v = 9),   file.path(tmp, "sub/solo.tsv"))

writeLines(c(
  "key\ttype\tlocal_subpath\tgs_uri",
  "multi\teQTL-pQTL\tsub/a.tsv\t-",
  "multi\tsQTL-pQTL\tsub/b.tsv\t-",
  "solo\t-\tsub/solo.tsv\t-",
  "nofile\t-\tsub/missing.tsv\t-"
), file.path(tmp, "inputs.tsv"))

MANIFEST <- file.path(tmp, "inputs.tsv")
INPUT    <- tmp
source("common/manifest.R")

# A manifest column named `key` must not shadow the `key` argument. This is the
# data.table non-standard-evaluation bug found during the Fig1 work: on a
# data.table, `m[m$key == key, ]` compares the column to itself and matches
# every row.
ok("multi resolves to exactly 2 rows", nrow(input_files("multi")) == 2L)
ok("solo resolves to exactly 1 row",   nrow(input_files("solo")) == 1L)
ok("input_path returns the solo path",
   input_path("solo") == file.path(tmp, "sub/solo.tsv"))

d <- load_input("multi")
ok("load_input row-binds every file for the key", nrow(d) == 5L)
ok("metadata attaches under its own manifest column name", "type" %in% names(d))
ok("metadata values follow their own file",
   identical(d$type, c("eQTL-pQTL", "eQTL-pQTL", "sQTL-pQTL", "sQTL-pQTL", "sQTL-pQTL")))
ok("reserved columns are never attached",
   !any(c("key", "local_subpath", "gs_uri", "path") %in% names(d)))

ok("'-' metadata is read as NA and not attached",
   !"type" %in% names(load_input("solo")))

expect_error("unknown key errors",                  input_files("nope"))
expect_error("multi-file key errors in input_path", input_path("multi"))
expect_error("missing file with no gs_uri errors",  input_files("nofile"))

cat("test_manifest.R: all assertions passed\n")
