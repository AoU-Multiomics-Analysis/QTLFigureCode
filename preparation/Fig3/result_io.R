# Local-only adapter for the existing figure result-table manifests.
# Unlike common/manifest.R, this adapter does not download or stage any input.
# Each run rebuilds scientific caches from result tables. Memoization is only
# within this R process; no old on-disk cache can bypass a preparation step.
prep_result_runtime <- function(figure, envir = parent.frame()) {
  manifest <- readr::read_tsv(file.path(figure, "inputs.tsv"), col_types = cols(.default = col_character()),
                             na = c("", "-"), show_col_types = FALSE)
  cache <- new.env(parent = emptyenv())
  input_files <- function(key) {
    rows <- manifest[manifest$key == key, , drop = FALSE]
    if (!nrow(rows)) stop("Unknown result-table key: ", key, call. = FALSE)
    rows$path <- map_chr(rows$local_subpath, function(relative) {
      if (key == "clingen") {
        root <- .prep_root("input")
        files <- list.files(file.path(root, relative), pattern = "[.]parquet$", full.names = FALSE)
        if (!length(files)) stop("No ClinGen parquet input under ", relative, call. = FALSE)
        return(prep_input("input", file.path(relative, sort(files)[1])))
      }
      candidates <- unique(c(relative, sub("splicing_SusieMerged[.]annotated", "splicing_SusieMerged..annotated", relative)))
      prep_input_any("input", candidates)
    })
    rows
  }
  input_path <- function(key) {
    rows <- input_files(key)
    if (nrow(rows) != 1L) stop("Expected one input for ", key, call. = FALSE)
    rows$path
  }
  load_input <- function(key, reader = data.table::fread, ..., attach_meta = TRUE) {
    rows <- input_files(key)
    metadata <- setdiff(names(rows), c("key", "local_subpath", "gs_uri", "source_repo", "source_workflow", "path"))
    map_dfr(seq_len(nrow(rows)), function(i) {
      data <- as_tibble(reader(rows$path[i], ...))
      if (attach_meta) for (column in metadata) if (!is.na(rows[[column]][i])) data[[column]] <- rows[[column]][i]
      data
    })
  }
  build_or_cache <- function(name, expr, refresh = FALSE) {
    if (!refresh && exists(name, envir = cache, inherits = FALSE)) return(cache[[name]])
    message("Prepare ", figure, "/", name, " from source results")
    value <- force(expr)
    cache[[name]] <- value
    prep_save_rds(value, "source", file.path(figure, "data", paste0(name, ".rds")))
    value
  }
  for (name in c("input_files", "input_path", "load_input", "build_or_cache")) assign(name, get(name), envir = envir)
  invisible(cache)
}
