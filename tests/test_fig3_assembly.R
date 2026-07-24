# Regression test for Figure 3's panel-tag lettering (see fig3-cleanup Finding 1
# and Finding 2). Panel tagging is a property of the ASSEMBLED figure, not of
# any single panel's code, so this builds the real panels/assembly logic
# against the fixture and inspects the actual tag grobs drawn by patchwork --
# reading assemble.R/supp_coloc_qc.R would not have caught the original bug,
# since patchwork's auto-tagging behavior (including the "inset counts as a
# tag" default, and the "ignore_tag on an inset also un-tags its own host
# plot when nested among siblings" edge case) is not visible from the source.
#
# Run from the project root:  Rscript tests/test_fig3_assembly.R
source("tests/fixtures.R")
fx <- make_fig3_fixture()
use_fig3_fixture(fx)

source("fig3/00_setup.R")
source("tests/helpers.R")

# Extract the tag letters actually drawn by a built patchwork figure.
#
# This is a two-step check, not a single heuristic:
#  1. Actually PRINT the figure to a null device and grid::grid.grab() the
#     resulting display list. This walks the exact same grob tree that would
#     hit a real graphics device (unlike inspecting the ggplot object graph,
#     which is how the original tagging bug survived code review), and grid
#     emits grobs in physical paint order, which for these layouts is reading
#     order (confirmed by rendering fig3_full.png / supp_coloc_qc while
#     developing this test).
#  2. Within that display list, a plot's tag is rendered as a bare, bold,
#     single-uppercase-letter text grob (`plot.tag`), and nothing else in
#     Figure 3 or the QC supplemental renders a single bare A-Z character
#     (every other text grob is an axis label, title, legend entry, or a
#     multi-character annotation) -- verified by dumping every text label in
#     the built figure while developing this test. So filtering the grabbed
#     display list down to single-character entries in LETTERS recovers
#     exactly the tag sequence, in draw order.
#
# Using patchwork:::patchworkGrob() (matching each plot's own gtable "tag"
# layout cell) was tried first and rejected: it reliably finds the right SET
# of tag grobs, but each nested sub-gtable uses its own local (t, l)
# coordinate frame, so sorting by (t, l) across nesting levels does not
# recover a single consistent reading order -- confirmed empirically.
extract_tag_labels <- function(patchwork_obj) {
  grDevices::pdf(NULL)
  on.exit({ grDevices::dev.off() })
  grid::grid.newpage()
  print(patchwork_obj)
  grabbed <- grid::grid.grab()

  labels <- character(0)
  walk <- function(x) {
    if (is.null(x)) return(invisible())
    if (inherits(x, "text")) {
      lab <- paste(x$label, collapse = "")
      if (nchar(lab) == 1L && lab %in% LETTERS) labels <<- c(labels, lab)
    } else if (inherits(x, "gtable")) {
      # gtable objects keep their children in $grobs, not $children.
      for (gr in x$grobs) walk(gr)
    } else if (!is.null(x$children)) {
      for (ch in x$children) walk(ch)
    }
  }
  walk(grabbed)
  labels
}

# --- Build the five panels exactly as fig3/run_all.R does ------------------
for (p in LETTERS[1:5]) source(sprintf("fig3/panels/panel%s.R", p))

P <- setNames(lapply(LETTERS[1:5],
                     function(x) readRDS(file.path(OBJ_DIR, paste0(x, ".rds")))),
              LETTERS[1:5])

coloc_counts <- P$A + P$B + plot_layout(widths = c(1, 1))
top_row      <- (coloc_counts | P$C) + plot_layout(widths = c(2.5, 1))
bottom_row   <- (P$D | P$E) + plot_layout(widths = c(1.3, 1))

Fig3 <- (top_row / bottom_row) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = SZ_TAG, face = "bold"))

fig3_tags <- extract_tag_labels(Fig3)

ok("Figure 3 carries exactly 5 tags", length(fig3_tags) == 5L)
ok("Figure 3's tags are exactly A, B, C, D, E in order",
   identical(fig3_tags, c("A", "B", "C", "D", "E")))
ok("no sixth tag (the inset) leaked into the sequence",
   !"F" %in% fig3_tags)

# --- Supplemental: exactly two tags, A on the density row, B on discovery --
source("fig3/supp_coloc_qc.R")

supp_tags <- extract_tag_labels(SuppColocQC)

ok("the QC supplemental carries exactly 2 tags", length(supp_tags) == 2L)
ok("the QC supplemental's tags are exactly A, B in order",
   identical(supp_tags, c("A", "B")))

cat("test_fig3_assembly.R: all assertions passed\n")
