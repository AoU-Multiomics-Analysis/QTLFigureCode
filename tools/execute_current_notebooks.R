#!/usr/bin/env Rscript
# Execute R notebook cells without a Python Jupyter installation.
suppressPackageStartupMessages({
  library(jsonlite)
  library(evaluate)
})

arguments <- commandArgs(trailingOnly = TRUE)
index <- read_json("notebooks/current/index.json", simplifyVector = FALSE)
targets <- if (length(arguments)) arguments else names(index)
stopifnot(all(targets %in% names(index)))
reports <- list()

for (key in targets) {
  message("Execute notebook: ", key)
  path <- index[[key]]$path
  nb <- read_json(path, simplifyVector = FALSE)
  env <- new.env(parent = globalenv())
  count <- 0L
  image_count <- 0L
  for (i in seq_along(nb$cells)) {
    if (nb$cells[[i]]$cell_type != "code") next
    count <- count + 1L
    message("  Cell ", count)
    outputs <- list()
    append_output <- function(value) outputs[[length(outputs) + 1L]] <<- value
    image_output <- function(file) {
      image_count <<- image_count + 1L
      append_output(list(output_type = "display_data", metadata = structure(list(), names = character()),
                         data = list("image/png" = base64_enc(readBin(file, "raw", file.info(file)$size)))))
      invisible(NULL)
    }
    options(qtl.notebook.display = image_output)
    plot_output <- function(plot) {
      temp <- tempfile(fileext = ".png")
      png(temp, width = 1500, height = 1100, res = 150, bg = "white")
      replayPlot(plot)
      dev.off()
      image_output(temp)
      unlink(temp)
      invisible(NULL)
    }
    text_output <- function(text) {
      if (nzchar(text)) append_output(list(output_type = "stream", name = "stdout", text = text))
      invisible(NULL)
    }
    handler <- new_output_handler(
      text = text_output, graphics = plot_output,
      message = function(x) text_output(conditionMessage(x)),
      warning = function(x) text_output(paste0("Warning: ", conditionMessage(x))),
      value = function(value, visible, envir) {
        if (visible && !is.null(value)) {
          if (inherits(value, c("ggplot", "patchwork", "grob", "gTree"))) print(value)
          else stop("A cell returned a raw value. Display a plot or a status message instead.")
        }
        invisible(NULL)
      })
    # Cairo measures Unicode text and the saved Arial labels correctly.
    # The default PDF device can replace characters before the PNG replay.
    measurement <- tempfile(fileext = ".pdf")
    cairo_pdf(measurement, width = 10, height = 7)
    dev.control(displaylist = "enable")
    tryCatch(
      evaluate(paste0(unlist(nb$cells[[i]]$source), collapse = ""), envir = env,
               stop_on_error = 2L, output_handler = handler, new_device = FALSE),
      finally = {dev.off(); unlink(measurement)})
    nb$cells[[i]]$outputs <- outputs
    nb$cells[[i]]$execution_count <- count
  }
  records <- as.list(env$input_records)
  for (record in records) {
    root <- if (record$root == "QTL_SOURCE_ROOT") env$source_root else env$paper_root
    file <- file.path(root, record$path)
    stopifnot(identical(unname(tools::md5sum(file)), record$md5))
  }
  stopifnot(image_count > 0L)
  nb$metadata$qtl_execution <- list(reference_date = "2026-11-02", R = as.character(getRversion()),
                                     cells = count, images = image_count,
                                     source_inputs_unchanged = TRUE,
                                     input_records = unname(records))
  for (i in seq_along(nb$cells)) {
    if (!length(nb$cells[[i]]$metadata)) {
      nb$cells[[i]]$metadata <- structure(list(), names = character())
    }
  }
  write_json(nb, path, auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
  reports[[key]] <- list(cells = count, images = image_count, source_files = length(records),
                         source_inputs_unchanged = TRUE)
  message("  Passed: ", count, " code cells; ", image_count, " images.")
}
options(qtl.notebook.display = NULL)
report_dir <- Sys.getenv("QTL_NOTEBOOK_OUTPUT", "notebook_outputs")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
write_json(reports, file.path(report_dir, "execution_report.json"),
           auto_unbox = TRUE, pretty = TRUE)
