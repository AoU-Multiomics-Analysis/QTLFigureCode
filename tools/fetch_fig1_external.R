#!/usr/bin/env Rscript
# =============================================================================
# Stage Figure 1's external comparison QTL data.
#
# Fig1 panel D counts eGenes/sGenes that are NOVEL relative to three external
# studies. Those five inputs are the only ones the pipeline cannot stage by
# itself, for two different reasons:
#
#   gtex_egenes / gtex_sgenes   public, but only distributed INSIDE a bundled
#                               archive, so the manifest cannot express them
#                               as a plain gs_uri -- this script extracts them
#   interval_eqtl / interval_sqtl   already carry a gs:// URI in the manifest
#   topmed_eqtl                     "     "        "        "
#
# The INTERVAL/TOPMed files would be staged automatically on first use by
# input_files(); they are included here so one command leaves panel D's whole
# external layer on disk.
#
# Usage (from the project root):
#   Rscript tools/fetch_fig1_external.R              # fetch what is missing
#   Rscript tools/fetch_fig1_external.R --force      # re-fetch everything
#   Rscript tools/fetch_fig1_external.R --gtex-only  # skip the gs:// inputs
#   FIG1_INPUT=/path/to/staging Rscript tools/fetch_fig1_external.R
#
# BANDWIDTH WARNING. GTEx v10 publishes per-tissue eGenes/sGenes only inside
# GTEx_Analysis_v10_{eQTL,sQTL}.tar (2.56 GB and 1.96 GB, and gzip-compressed
# despite the .tar name). This script STREAMS each archive and writes only the
# one member it needs, so nothing large is stored -- but up to ~4.5 GB may be
# transferred. Members are NOT stored in alphabetical order, so tar's early exit
# saves an unpredictable amount: anywhere from most of the archive to none of it.
# Budget for the worst case on a metered connection.
#
# The archives' internal directory is GTEx_Analysis_v10_{e,s}QTL_updated/ -- not
# the name of the .tar -- which is why members are matched by GLOB below rather
# than by a hard-coded path.
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })

args     <- commandArgs(trailingOnly = TRUE)
force    <- "--force" %in% args
gtex_only<- "--gtex-only" %in% args

INPUT    <- Sys.getenv("FIG1_INPUT", ".")
MANIFEST <- Sys.getenv("FIG1_MANIFEST", "Fig1/inputs.tsv")
if (!file.exists(MANIFEST)) stop("manifest not found: ", MANIFEST)

# as.data.frame() defuses data.table's non-standard evaluation: the manifest has
# a column literally named `key`, so inside man[...] a `key` ARGUMENT resolves to
# the column instead, comparing it to itself and matching every row. This is the
# same trap common/manifest.R documents; it silently returns all 47 rows.
man <- as.data.frame(fread(MANIFEST, sep = "\t", na.strings = c("", "-"),
                           colClasses = "character"))
row_for <- function(key) {
  r <- man[man$key == key, ]
  if (nrow(r) != 1L) stop("expected exactly 1 manifest row for '", key, "', got ", nrow(r))
  r
}

GTEX_BASE <- paste0("https://storage.googleapis.com/adult-gtex/bulk-qtl/v10/",
                    "single-tissue-cis-qtl/")
# key -> (archive, the member to pull out). Matched as a GLOB so the archive's
# internal directory name does not have to be hard-coded.
GTEX <- list(
  gtex_egenes = list(tar = "GTEx_Analysis_v10_eQTL.tar",
                     member = "*Whole_Blood.v10.eGenes.txt.gz"),
  gtex_sgenes = list(tar = "GTEx_Analysis_v10_sQTL.tar",
                     member = "*Whole_Blood.v10.sGenes.txt.gz"))

# GNU tar and bsdtar spell "stop once the pattern is matched" differently, and
# bsdtar needs no --wildcards. Getting this wrong silently extracts the whole
# archive to disk, which is the thing this script exists to avoid.
tar_flags <- function() {
  is_gnu <- any(grepl("GNU tar", suppressWarnings(system2("tar", "--version",
                       stdout = TRUE, stderr = TRUE)), fixed = TRUE))
  if (is_gnu) c("--wildcards", "--occurrence=1") else c("--fast-read")
}

ok_size <- function(p) file.exists(p) && file.info(p)$size > 0

fetch_gtex <- function(key) {
  r    <- row_for(key)
  dest <- file.path(INPUT, r$local_subpath)
  if (!force && ok_size(dest)) { message("  cached  ", r$local_subpath); return(TRUE) }
  spec <- GTEX[[key]]
  url  <- paste0(GTEX_BASE, spec$tar)
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)

  tmp <- file.path(tempdir(), paste0("gtex_", key))
  unlink(tmp, recursive = TRUE); dir.create(tmp, recursive = TRUE)

  message("  fetching ", spec$tar, " (streamed; extracting ", basename(spec$member), ")")
  # curl is killed by SIGPIPE as soon as tar has what it wants, so the pipeline's
  # exit status is not meaningful -- verify by looking for the file instead.
  cmd <- sprintf("curl -sL --fail %s | tar -xzf - -C %s %s %s",
                 shQuote(url), shQuote(tmp),
                 paste(tar_flags(), collapse = " "), shQuote(spec$member))
  suppressWarnings(system(cmd, ignore.stderr = TRUE))

  hit <- list.files(tmp, pattern = glob2rx(basename(spec$member)),
                    recursive = TRUE, full.names = TRUE)
  if (!length(hit)) {
    message("  FAILED  ", r$local_subpath,
            " -- '", spec$member, "' not found in ", spec$tar)
    return(FALSE)
  }
  file.copy(hit[1], dest, overwrite = TRUE)
  unlink(tmp, recursive = TRUE)
  message("  wrote   ", r$local_subpath, "  (",
          format(file.info(dest)$size / 1e6, digits = 3), " MB)")
  TRUE
}

fetch_gs <- function(key) {
  r    <- row_for(key)
  dest <- file.path(INPUT, r$local_subpath)
  if (!force && ok_size(dest)) { message("  cached  ", r$local_subpath); return(TRUE) }
  if (is.na(r$gs_uri)) {
    message("  SKIP    ", r$local_subpath, " -- no gs_uri in the manifest")
    return(FALSE)
  }
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  message("  fetching ", r$gs_uri)
  st <- suppressWarnings(system2("gsutil", c("cp", r$gs_uri, dest)))
  if (st != 0L || !ok_size(dest)) {
    message("  FAILED  ", r$local_subpath,
            " -- gsutil cp returned ", st,
            ". This bucket is controlled-access; check `gcloud auth login`.")
    return(FALSE)
  }
  message("  wrote   ", r$local_subpath, "  (",
          format(file.info(dest)$size / 1e6, digits = 3), " MB)")
  TRUE
}

message("staging root: ", normalizePath(INPUT, mustWork = FALSE))
message("\nGTEx v10 (public, extracted from the bundled archive):")
res <- vapply(names(GTEX), fetch_gtex, logical(1))

if (!gtex_only) {
  message("\nINTERVAL / TOPMed (gs:// from the manifest):")
  res <- c(res, vapply(c("interval_eqtl", "interval_sqtl", "topmed_eqtl"),
                       fetch_gs, logical(1)))
}

message("\n", sum(res), "/", length(res), " staged")
if (!all(res)) {
  message("Still missing:")
  for (k in names(res)[!res])
    message("  ", k, " -> ", file.path(INPUT, row_for(k)$local_subpath))
  quit(save = "no", status = 1)
}
