# Asserts the assembled Figure 4 carries exactly the four tags A, B, C, D, by
# rendering and inspecting the tag grobs (not by reading the code) -- the same
# check that caught Fig3's mislettered inset. Requires MatchIt (Panel C) and
# arrow (Panel D). Run from the project root.
source("tests/fixtures.R")
source("tests/helpers.R")
if (!requireNamespace("MatchIt", quietly = TRUE) ||
    !requireNamespace("arrow", quietly = TRUE)) {
  cat("  (skipped test_fig4_assembly: MatchIt/arrow not installed)\n")
  quit(status = 0)
}
fx <- make_fig4_fixture(); use_fig4_fixture(fx)

# Proven extractor, copied from tests/test_fig3_assembly.R: render the built
# figure to a null device, grab the grid display list, and collect the bare
# single-uppercase-letter text grobs (the plot.tag grobs). Nothing else in the
# figure draws a bare A-Z character.
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

for (p in LETTERS[1:4]) source(sprintf("fig4/panels/panel%s.R", p))
source("fig4/panels/panelA_labeled.R")
source("fig4/assemble.R")

tags <- extract_tag_labels(Fig4)
# The layout wraps the B/C row in wrap_elements, so draw order is not
# necessarily A,B,C,D; assert the SET is exactly the four tags, four in total,
# with none missing, extra, or duplicated.
ok("Figure 4 carries exactly 4 tags", length(tags) == 4L)
ok("Figure 4's tags are exactly A, B, C, D",
   identical(sort(unique(tags)), c("A", "B", "C", "D")))

# The labelled-A variant shares compose_fig4(), so it must tag identically. The
# shorthand labels are multi-character and so cannot be mistaken for tag grobs.
labelled_tags <- extract_tag_labels(Fig4_labeled)
ok("Figure 4 (labelled A) carries exactly 4 tags", length(labelled_tags) == 4L)
ok("Figure 4 (labelled A) tags are exactly A, B, C, D",
   identical(sort(unique(labelled_tags)), c("A", "B", "C", "D")))

cat("test_fig4_assembly.R: all assertions passed\n")
