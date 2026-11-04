# Check portable paths and clear failures before any figure input is read.
previous <- Sys.getenv(c("QTL_PAPER_ROOT", "QTL_SOURCE_ROOT", "QTL_NOTEBOOK_OUTPUT", "QTL_PREP_ROOT"))
root <- tempfile("qtl input with spaces ")
dir.create(root)
writeLines("fixture", file.path(root, "input.txt"))
Sys.setenv(QTL_PAPER_ROOT = root, QTL_SOURCE_ROOT = root,
           QTL_NOTEBOOK_OUTPUT = file.path(root, "outputs"))
env <- new.env(parent = globalenv())
sys.source("notebooks/current/00_setup.R", envir = env)
before <- tools::md5sum(file.path(root, "input.txt"))
stopifnot(identical(env$paper_file("input.txt"), normalizePath(file.path(root, "input.txt"))))
stopifnot(identical(env$source_file("input.txt"), normalizePath(file.path(root, "input.txt"))))
# Prepared output files must not replace the original notebook inputs.
prepared <- file.path(root, "prepared")
dir.create(file.path(prepared, "paper"), recursive = TRUE)
writeLines("prepared fixture", file.path(prepared, "paper", "input.txt"))
Sys.setenv(QTL_PREP_ROOT = prepared)
overlay <- new.env(parent = globalenv())
sys.source("notebooks/current/00_setup.R", envir = overlay)
stopifnot(identical(overlay$paper_file("input.txt"), normalizePath(file.path(root, "input.txt"))))
stopifnot(identical(overlay$source_file("input.txt"), normalizePath(file.path(root, "input.txt"))))
stopifnot(all(vapply(as.list(overlay$input_records), function(x) x$root != "QTL_PREP_ROOT", logical(1))))
message_for <- function(expression) tryCatch({force(expression); ""}, error = conditionMessage)
stopifnot(grepl("Missing or unreadable input", message_for(env$paper_file("absent.tsv")), fixed = TRUE))
env$paper_root <- "gs://bucket/input"
stopifnot(grepl("Input localization error", message_for(env$paper_file("x")), fixed = TRUE))
env$paper_root <- ""
stopifnot(grepl("Set QTL_PAPER_ROOT", message_for(env$paper_file("x")), fixed = TRUE))
target <- env$output_file("nested/result.txt")
stopifnot(dir.exists(dirname(target)), !file.exists(target),
          identical(before, tools::md5sum(file.path(root, "input.txt"))),
          length(as.list(env$input_records)) == 2L)
unlink(root, recursive = TRUE)
for (name in names(previous)) {
  if (nzchar(previous[[name]])) do.call(Sys.setenv, setNames(list(previous[[name]]), name))
  else Sys.unsetenv(name)
}
message("Current notebook input tests passed.")
