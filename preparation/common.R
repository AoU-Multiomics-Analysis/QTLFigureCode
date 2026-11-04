# Local input roots and isolated outputs for figure-data preparation.
# Run preparation scripts from the repository root.
suppressPackageStartupMessages(library(tidyverse))
.prep_state <- new.env(parent = emptyenv())
.prep_state$inputs <- list()
.prep_state$outputs <- list()

prep_init <- function() {
  source_root <- Sys.getenv("QTL_SOURCE_ROOT", getwd())
  input_root <- Sys.getenv("QTL_INPUT_ROOT", source_root)
  config <- list(paper_root = Sys.getenv("QTL_PAPER_ROOT"), source_root = source_root,
                 input_root = input_root,
                 susie_root = Sys.getenv("QTL_SUSIE_ROOT", input_root),
                 output_root = Sys.getenv("QTL_PREP_ROOT", file.path(getwd(), "notebook_outputs", "prepared")))
  if (grepl("^[[:alpha:]][[:alnum:]+.-]*://", config$output_root)) {
    stop("The preparation output root must be a local directory.", call. = FALSE)
  }
  output <- normalizePath(config$output_root, mustWork = FALSE)
  for (root in config[c("paper_root", "source_root", "input_root", "susie_root")]) {
    if (nzchar(root) && identical(output, normalizePath(root, mustWork = FALSE))) {
      stop("The preparation output root must differ from every input root.", call. = FALSE)
    }
  }
  .prep_state$config <- config
  config
}

.prep_relative <- function(relative) {
  if (length(relative) != 1L || is.na(relative) || !nzchar(relative) ||
      grepl("^/|^[A-Za-z]:|://|(^|[/\\\\])\\.\\.([/\\\\]|$)", relative)) {
    stop("Use a nonempty relative file path without parent traversal.", call. = FALSE)
  }
  relative
}

.prep_root <- function(name) {
  key <- c(paper = "paper_root", source = "source_root", input = "input_root", susie = "susie_root")[[name]]
  if (is.null(key)) stop("Unknown input root: ", name, call. = FALSE)
  value <- .prep_state$config[[key]]
  if (!nzchar(value)) stop("Set the local ", name, " input root before preparation.", call. = FALSE)
  if (grepl("^[[:alpha:]][[:alnum:]+.-]*://", value)) {
    stop("Input localization error: ", name, " must be a local directory.", call. = FALSE)
  }
  value
}

.prep_remember <- function(path, root, relative) {
  if (!file.exists(path) || dir.exists(path) || file.access(path, 4) != 0) {
    stop("Missing or unreadable input: ", relative, " (", root, ").", call. = FALSE)
  }
  path <- normalizePath(path, mustWork = TRUE)
  if (is.null(.prep_state$inputs[[path]])) {
    .prep_state$inputs[[path]] <- list(root = root, path = relative,
                                     md5 = unname(tools::md5sum(path)))
  }
  path
}

prep_input <- function(root_name, relative) {
  relative <- .prep_relative(relative)
  if (root_name %in% c("paper", "source")) {
    generated <- file.path(.prep_state$config$output_root, root_name, relative)
    if (file.exists(generated)) return(.prep_remember(generated, "prepared", file.path(root_name, relative)))
  }
  .prep_remember(file.path(.prep_root(root_name), relative), root_name, relative)
}

# Use original-root files only when checking against the accepted result.
prep_reference <- function(root_name, relative) {
  relative <- .prep_relative(relative)
  .prep_remember(file.path(.prep_root(root_name), relative), root_name, relative)
}

prep_has_reference <- function(root_name, relative) {
  relative <- .prep_relative(relative)
  root <- tryCatch(.prep_root(root_name), error = function(e) NULL)
  if (is.null(root)) return(FALSE)
  path <- file.path(root, relative)
  file.exists(path) && !dir.exists(path) && file.access(path, 4) == 0
}

prep_root <- function(root_name) .prep_root(root_name)

prep_input_any <- function(root_name, relative_candidates) {
  for (relative in relative_candidates) {
    .prep_relative(relative)
    generated <- file.path(.prep_state$config$output_root, root_name, relative)
    source <- file.path(.prep_root(root_name), relative)
    if ((root_name %in% c("paper", "source") && file.exists(generated)) || file.exists(source)) {
      return(prep_input(root_name, relative))
    }
  }
  stop("Missing input: none of ", paste(relative_candidates, collapse = ", "),
       " is readable under the ", root_name, " root.", call. = FALSE)
}

prep_input_file <- function(variable, fallback_root = "input", relative) {
  override <- Sys.getenv(variable)
  if (!nzchar(override)) return(prep_input_any(fallback_root, relative))
  if (grepl("^[[:alpha:]][[:alnum:]+.-]*://", override)) {
    stop("Input localization error: ", variable, " must name a local file.", call. = FALSE)
  }
  .prep_remember(override, paste0("env:", variable), basename(override))
}

prep_output <- function(namespace, relative) {
  if (!namespace %in% c("paper", "source")) {
    stop("Output namespace must be paper or source.", call. = FALSE)
  }
  relative <- .prep_relative(relative)
  path <- file.path(.prep_state$config$output_root, namespace, relative)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  path <- normalizePath(path, mustWork = FALSE)
  if (!is.null(.prep_state$inputs[[path]])) {
    stop("Output would overwrite an input: ", relative, call. = FALSE)
  }
  .prep_state$outputs[[path]] <- list(root = namespace, path = relative)
  path
}

prep_write_tsv <- function(data, namespace, relative) {
  readr::write_tsv(data, prep_output(namespace, relative), quote = "needed")
  message("Wrote derived table: ", namespace, "/", relative)
  invisible(data)
}

prep_save_rds <- function(object, namespace, relative) {
  saveRDS(object, prep_output(namespace, relative))
  message("Wrote derived object: ", namespace, "/", relative)
  invisible(object)
}

prep_finish <- function(key) {
  for (path in names(.prep_state$inputs)) {
    if (!identical(unname(tools::md5sum(path)), .prep_state$inputs[[path]]$md5)) {
      stop("An input changed during preparation: ", .prep_state$inputs[[path]]$path, call. = FALSE)
    }
  }
  folder <- file.path(.prep_state$config$output_root, "logs")
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  report <- list(step = key, source_inputs_unchanged = TRUE,
                 inputs = unname(.prep_state$inputs), outputs = unname(.prep_state$outputs),
                 R = as.character(getRversion()))
  jsonlite::write_json(report, file.path(folder, paste0(key, ".json")),
                       auto_unbox = TRUE, pretty = TRUE, null = "null")
  message("Preparation passed: ", key, "; source inputs unchanged.")
  invisible(report)
}

invisible(prep_init())
