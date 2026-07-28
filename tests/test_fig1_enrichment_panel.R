# Fig1 panel H (rendered tag I) — functional-enrichment panel rendering.
# Run from the project root:  Rscript tests/test_fig1_enrichment_panel.R
#
# Fisher's exact test returns an infinite odds ratio whenever a cell of the 2x2
# table is empty, which really happens for the sparsest pQTL annotations. Those
# points used to be handed to ggplot as +Inf, which draws them on the panel edge
# clipped in half. RenderFinemappingEnrichment() is the pure seam that has to
# keep every marker inside the panel without inventing a finite estimate.
source("tests/helpers.R")
suppressMessages(source("fig1/00_setup.R"))

# Two functional groups; the pQTL "stop gained" row is the unbounded one.
enrichment <- data.frame(
  name             = c("missense_variant", "missense_variant", "stop_gained",
                       "PLS", "dELS"),
  type             = c("eQTL", "pQTL", "pQTL", "eQTL", "eQTL"),
  odds_ratio       = c(10, 47, Inf, 4, 2),
  conf.low         = c(5, 33, 15, 3, 1.5),
  conf.high        = c(20, 68, Inf, 6, 3),
  functional_group = c("Coding Annotation", "Coding Annotation",
                       "Coding Annotation", "Regulatory Annotation",
                       "Regulatory Annotation"),
  annotation_clean = c("Missense", "Missense", "Stop gained",
                       "Promoter", "Distal enhancer"),
  row              = c(1, 1, 2, 2, 1),
  stringsAsFactors = FALSE
)

p  <- RenderFinemappingEnrichment(enrichment)
pd <- ggplot_build(p)

# Layer 3 is the size-7 geom_point; layer 2 is the geom_pointrange.
points <- pd$data[[3]]
ranges <- pd$data[[2]]

ok("no plotted point is infinite", all(is.finite(points$y)))
ok("no confidence interval runs to infinity",
   all(is.finite(ranges$ymin)) && all(is.finite(ranges$ymax)))

# The panel is log10-transformed, so build-stage values are already log10.
cap <- max(c(10, 47, 4, 2, 20, 68, 6, 3)) * 1.5    # 1.5 x the largest bounded value
ok("the unbounded estimate is drawn at the ceiling",
   isTRUE(all.equal(max(points$y), log10(cap))))
ok("bounded estimates keep their own values",
   isTRUE(all.equal(sort(points$y[points$y < log10(cap)]),
                    sort(log10(c(10, 47, 4, 2))))))
ok("the unbounded interval is capped at the same ceiling",
   isTRUE(all.equal(max(ranges$ymax), log10(cap))))

ok("the unbounded estimate is drawn as a triangle, not a circle",
   points$shape[which.max(points$y)] == 24)
ok("every bounded estimate stays a circle",
   all(points$shape[points$y < log10(cap)] == 21))

# Clipping is the actual defect being guarded against: a marker exactly on the
# panel edge is drawn as a half circle.
y_range <- pd$layout$panel_params[[1]]$y.range
ok("every marker sits strictly inside the panel, not on its edge",
   all(points$y > y_range[1]) && all(points$y < y_range[2]))
ok("the no-enrichment reference line is inside the axis range",
   y_range[1] < 0 && y_range[2] > 0)   # log10(1) == 0

# The grey striping is drawn with ymax = Inf so it fills the panel; that must not
# feed back into the axis range. Headroom above the topmost marker should be a
# few percent of the panel -- enough to clear a size-7 marker, no more.
headroom <- (y_range[2] - max(points$y)) / diff(y_range)
ok("there is headroom above the topmost marker", headroom > 0.01)
ok("the shading bands do not stretch the axis", headroom < 0.15)

# --- the panel as it is actually built: PIP logistic regression ---------------
# Panel H now comes from RenderGlmEnrichment(). glm() returns a log-odds slope
# and the panel shows it exponentiated on a log10 axis, so nothing here can be
# infinite -- but the log transform brings its own trap, below.
glm_results <- data.frame(
  name             = c("missense_variant", "missense_variant", "stop_gained",
                       "PLS", "dELS"),
  type             = c("eQTL", "pQTL", "pQTL", "eQTL", "eQTL"),
  estimate         = c(1.34, 3.03, 4.80, 2.25, 0.88),
  std.error        = c(0.03, 0.07, 0.62, 0.06, 0.02),
  conf.low         = c(1.28, 2.89, 3.58, 2.13, 0.84),
  conf.high        = c(1.40, 3.17, 6.02, 2.37, 0.92),
  padj             = c(1e-40, 1e-60, 1e-14, 1e-90, 1e-30),
  functional_group = c("Coding Annotation", "Coding Annotation",
                       "Coding Annotation", "Regulatory Annotation",
                       "Regulatory Annotation"),
  annotation_clean = c("Missense", "Missense", "Stop gained",
                       "Promoter", "Distal enhancer"),
  row              = c(1, 1, 2, 2, 1),
  stringsAsFactors = FALSE
)

gp  <- RenderGlmEnrichment(glm_results)
gpd <- ggplot_build(gp)

# Selected by geom rather than by position: adding or removing a layer silently
# renumbers positional indices, and these assertions then check the wrong data.
layer_data_for <- function(build, geom_class) {
  idx <- which(vapply(build$plot$layers,
                      function(l) inherits(l$geom, geom_class), logical(1)))
  if (length(idx) != 1L) stop("expected exactly one ", geom_class, " layer")
  build$data[[idx]]
}
g_ranges <- layer_data_for(gpd, "GeomLinerange")
g_points <- layer_data_for(gpd, "GeomPoint")
g_rules  <- layer_data_for(gpd, "GeomSegment")

ok("the slope is exponentiated onto the odds-ratio scale",
   isTRUE(all.equal(sort(g_points$y), sort(log10(exp(glm_results$estimate))))))
ok("confidence bounds are exponentiated with it",
   isTRUE(all.equal(sort(g_ranges$ymin), sort(log10(exp(glm_results$conf.low))))))
ok("every point is finite", all(is.finite(g_points$y)))

g_yrange <- gpd$layout$panel_params[[1]]$y.range
ok("the no-enrichment line at an odds ratio of 1 is inside the axis range",
   g_yrange[1] < 0 && g_yrange[2] > 0)   # log10(1) == 0
ok("every marker sits strictly inside the panel",
   all(g_points$y > g_yrange[1]) && all(g_points$y < g_yrange[2]))

# The grey rules replace the x axis line, so there must be exactly one per
# annotation and it must sit below the reference line rather than on it --
# taking their height from the data alone put them on top of it.
n_annotations <- nrow(dplyr::distinct(glm_results, functional_group, annotation_clean))
ok("one grey rule per annotation", nrow(g_rules) == n_annotations)
ok("the rules sit below the odds-ratio-of-1 line", all(g_rules$y < 0))
ok("the rules are inside the panel", all(g_rules$y > g_yrange[1]))
ok("the x axis line is suppressed in favour of the rules",
   inherits(gp$theme$axis.line.x, "element_blank"))

# --- assay colours -----------------------------------------------------------
# Panels D, G and H all colour by assay but list their categories in different
# orders (eGene/sGene/pGene vs Expression/Protein/Splicing vs eQTL/pQTL/sQTL).
# Taking colours positionally off the palette therefore gave the splicing bar
# the protein colour in one panel and not in the others.
ok("every alias for expression shares one colour",
   length(unique(AssayColors[c("Expression", "eQTL", "eGene")])) == 1)
ok("every alias for protein shares one colour",
   length(unique(AssayColors[c("Protein", "pQTL", "pGene")])) == 1)
ok("every alias for splicing shares one colour",
   length(unique(AssayColors[c("Splicing", "sQTL", "sGene")])) == 1)
ok("the three assays are given three different colours",
   length(unique(AssayColors[c("eGene", "pGene", "sGene")])) == 3)

# The enrichment panel is the reference for which colour means what.
fills <- pd$data[[3]][, c("fill")]
types <- c("eQTL", "pQTL", "pQTL", "eQTL", "eQTL")   # row order of `enrichment`
ok("the enrichment panel draws each assay in its assigned colour",
   all(tolower(fills) == tolower(AssayColors[types])))

cat("test_fig1_enrichment_panel.R: all assertions passed\n")
