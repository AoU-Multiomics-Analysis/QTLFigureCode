# Fig5 rendered-tag test. Patchwork's tagging of a heterogeneous ggplot /
# patchwork-with-inset / cowplot layout is only visible in the RENDERED output,
# so we build stand-in panels matching each real panel's structure, run the real
# assemble.R, and walk the drawn grob tree for the tag letters. Run from root.
source("tests/helpers.R")
suppressPackageStartupMessages(library(patchwork))

root <- file.path(tempdir(), "fig5_assembly")
unlink(root, recursive = TRUE)
dir.create(file.path(root, "objects"), recursive = TRUE, showWarnings = FALSE)
Sys.setenv(FIG5_OBJ = file.path(root, "objects"),
           FIG5_DATA = file.path(root, "data"),
           FIG5_OUTDIR = root)
source("fig5/00_setup.R")   # leaves FIG5_MANIFEST/INPUT default -> committed reference tables

# The proven extractor from tests/test_fig3_assembly.R: render to a null device,
# grid-grab the display list, and collect single-letter text grobs in draw order.
extract_tag_labels <- function(patchwork_obj) {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  grid::grid.newpage(); print(patchwork_obj)
  grabbed <- grid::grid.grab()
  labels <- character(0)
  walk <- function(x) {
    if (is.null(x)) return(invisible())
    if (inherits(x, "text")) {
      lab <- paste(x$label, collapse = "")
      if (nchar(lab) == 1L && lab %in% LETTERS) labels <<- c(labels, lab)
    } else if (inherits(x, "gtable")) { for (gr in x$grobs) walk(gr) }
    else if (!is.null(x$children)) { for (ch in x$children) walk(ch) }
  }
  walk(grabbed); labels
}

blank <- function() ggplot() + theme_void(base_family = FONT)
inset_pair <- function()   # main-with-inset over a gene-track stand-in (like D/F)
  (blank() + inset_element(blank(), 0.2, 0.2, 0.9, 0.9, ignore_tag = TRUE)) / blank()

standins <- list(
  A = blank() / blank(),          # 2-plot patchwork (ancestry / ClinVar)
  B = blank(), C = blank(), E = blank(),
  D = inset_pair(), F = inset_pair(),
  G = wrap_elements(full = blank()))
for (nm in names(standins))
  saveRDS(standins[[nm]], file.path(OBJ_DIR, paste0(nm, ".rds")))

source("fig5/assemble.R")

tags <- extract_tag_labels(Fig5)
ok("Figure 5 carries exactly 7 tags", length(tags) == 7L)
ok("its tags are exactly A, B, C, D, E, F, G",
   setequal(unique(tags), c("A", "B", "C", "D", "E", "F", "G")))
ok("Figure5.pdf was written", file.exists(file.path(root, "Figure5.pdf")))

cat("test_fig5_assembly.R passed.\n")
