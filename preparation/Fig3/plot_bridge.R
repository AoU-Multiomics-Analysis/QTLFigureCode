# Preserve accepted layouts while replacing their scientific data with values
# rebuilt by preparation steps. A failed comparison stops before any plot is
# saved. Templates provide layout, scales, and fixed graphical annotations.
suppressPackageStartupMessages(library(patchwork))
prepared_file <- function(namespace, relative) {
  .prep_relative(relative)
  path <- file.path(.prep_state$config$output_root, namespace, relative)
  .prep_remember(path, "prepared", file.path(namespace, relative))
}
prepared_rds <- function(namespace, relative) readRDS(prepared_file(namespace, relative))
prepared_tsv <- function(namespace, relative) readr::read_tsv(prepared_file(namespace, relative), show_col_types = FALSE)

replace_prepared_frame <- function(reference, data, label) {
  stopifnot(is.data.frame(reference), is.data.frame(data))
  absent <- setdiff(names(reference), names(data))
  if (length(absent)) stop(label, ": missing prepared columns: ", paste(absent, collapse = ", "), call. = FALSE)
  data <- data[, names(reference), drop = FALSE]
  for (column in names(reference)) {
    if (is.factor(reference[[column]])) data[[column]] <- factor(data[[column]], levels = levels(reference[[column]]))
  }
  if (!isTRUE(all.equal(reference, data, tolerance = 1e-10, check.attributes = FALSE))) {
    stop(label, ": prepared scientific values differ from the accepted layout reference; review the source inputs.", call. = FALSE)
  }
  message("Verified and replaced scientific data: ", label)
  data
}

replace_prepared_subset <- function(reference, data, keys, label) {
  index <- data %>% mutate(.prepared_row = row_number()) %>% select(all_of(keys), .prepared_row)
  match <- reference %>% select(all_of(keys)) %>% left_join(index, by = keys, relationship = "many-to-one")
  stopifnot(!anyNA(match$.prepared_row))
  replace_prepared_frame(reference, data[match$.prepared_row, , drop = FALSE], label)
}
