# =============================================================================
# Shared input manifest and lazy gsutil staging. Sourced by every figure's
# 00_setup.R, which must define these two globals FIRST:
#   MANIFEST  path to the figure's TSV manifest
#   INPUT     staging root (the cwd by default)
#
# Manifest columns `key`, `local_subpath` and `gs_uri` are reserved. EVERY OTHER
# COLUMN is metadata, attached by load_input() to the rows it returns, under its
# own name. Empty or "-" is read as NA and never attached.
#
# Files are cached at INPUT/<local_subpath>. A gs_uri is staged ONLY when the
# local copy is missing (gs:// via gsutil; http(s):// or file:// via
# download.file — see .stage_uri), so a pre-staged file is never re-fetched.
# =============================================================================
# Columns the loader treats as manifest STRUCTURE, never as data. Anything else
# is a metadata column and load_input() attaches it to every row of the loaded
# table -- so provenance columns must live here, or `source_repo` would appear
# as a column inside every SuSiE / enloc table Fig1-5 read.
.MANIFEST_RESERVED <- c("key", "local_subpath", "gs_uri",
                        "source_repo", "source_workflow")

.manifest_env <- new.env(parent = emptyenv())

input_manifest <- function() {
  # Re-read when MANIFEST changes so several figures (or tests) can share a session.
  if (is.null(.manifest_env$m) || !identical(.manifest_env$path, MANIFEST)) {
    # as_tibble() avoids data.table non-standard evaluation: the manifest has a
    # column literally named "key", which would shadow a `key` filter argument.
    .manifest_env$m <- tibble::as_tibble(
      data.table::fread(MANIFEST, sep = "\t", na.strings = c("", "-"),
                        colClasses = "character"))
    .manifest_env$path <- MANIFEST
  }
  .manifest_env$m
}

manifest_rows <- function(key) {
  m <- input_manifest()
  rows <- m[m$key == key, ]
  if (nrow(rows) == 0L) stop("Unknown input key '", key, "' in ", MANIFEST)
  rows
}

# Stage a single remote URI to a local destination, dispatching on scheme:
#   gs://            -> gsutil cp   (needs the Google Cloud SDK)
#   http(s):// file:// -> download.file
# Errors (naming the URI) if the scheme is unsupported or the fetch fails.
.stage_uri <- function(uri, dest) {
  if (grepl("^gs://", uri)) {
    status <- system2("gsutil", c("cp", uri, dest))
    if (status != 0L || !file.exists(dest)) stop("gsutil cp failed for ", uri)
  } else if (grepl("^(https?|file)://", uri)) {
    ok <- tryCatch({ utils::download.file(uri, dest, mode = "wb", quiet = TRUE); TRUE },
                   error = function(e) FALSE)
    if (!isTRUE(ok) || !file.exists(dest)) stop("download.file failed for ", uri)
  } else {
    stop("Unsupported URI scheme for staging (key's gs_uri = '", uri, "').")
  }
}

# Ensure every file for `key` exists locally (staging from gs_uri if needed) and
# return the manifest rows plus a resolved local `path` column.
input_files <- function(key) {
  rows  <- manifest_rows(key)
  paths <- file.path(INPUT, rows$local_subpath)
  for (i in seq_len(nrow(rows))) {
    lp <- paths[i]
    if (file.exists(lp)) next
    gs <- rows$gs_uri[i]
    if (is.na(gs))
      stop("Missing local input '", lp, "' and no gs_uri to stage it (key '", key, "').")
    dir.create(dirname(lp), showWarnings = FALSE, recursive = TRUE)
    message("  staging ", gs, " -> ", lp)
    .stage_uri(gs, lp)
  }
  rows$path <- paths
  tibble::as_tibble(rows)
}

# Single-file key -> its local path (staged if needed).
input_path <- function(key) {
  f <- input_files(key)
  if (nrow(f) != 1L)
    stop("Key '", key, "' maps to ", nrow(f), " files; use input_files()/load_input().")
  f$path
}

# Metadata columns of a resolved input_files() tibble.
manifest_meta_cols <- function(f) setdiff(names(f), c(.MANIFEST_RESERVED, "path"))

# Read every file for `key` and row-bind, attaching each metadata column under
# its own name. `reader` defaults to fread; override it (and pass extra args via
# ...) for other formats.
load_input <- function(key, reader = data.table::fread, ..., attach_meta = TRUE) {
  f    <- input_files(key)
  meta <- manifest_meta_cols(f)
  purrr::map_dfr(seq_len(nrow(f)), function(i) {
    d <- tibble::as_tibble(reader(f$path[i], ...))
    if (attach_meta) {
      for (mc in meta) {
        v <- f[[mc]][i]
        if (!is.na(v)) d[[mc]] <- v
      }
    }
    d
  })
}
