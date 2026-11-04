# Real fixture checks for plots rebuilt from source-derived scientific frames.
suppressPackageStartupMessages({library(tidyverse); library(patchwork)})
source("preparation/Fig3/methods.R")
source("preparation/Fig3/symmetric_methods.R")
source("preparation/Fig3/current_methods.R")
expect <- function(value, label) {
  if (!isTRUE(value)) stop("FAIL: ", label, call. = FALSE)
  message("PASS: ", label)
}
renderer <- "notebooks/current/plots/Fig3.R"
expect(file.exists(renderer), "Figure 3 has source-data plot builders")
source(renderer)

# A uses the shared-gene denominator, with three pair calls per gene.
gene <- crossing(Gene = c("shared1", "shared2", "excluded"),
                 type = c("eQTL-pQTL", "eQTL-sQTL", "sQTL-pQTL")) %>%
  mutate(GLCP = c(.99, .2, .2, .2, .2, .2, .99, .99, .2), Coloc = GLCP > .95)
a <- fig3_gene_counts(gene, c("shared1", "shared2"))
pa <- fig3_panel_a(a)
expect(sum(pa$data$n) == 2L && sum(pa$data$prop) == 1,
       "A retains the shared-gene denominator")
expect(isTRUE(all.equal(ggplot_build(pa)$data[[1]]$y, a$prop)),
       "A renders the recalculated gene proportions")

# B counts graph components per gene, but restricts only the inset denominator.
components <- tibble(Gene = c("shared1", "shared1", "other", "shared2"),
                     cs_id = c("a", "b", "c", "d"), type = c("eQTL", "pQTL", "sQTL", "eQTL"),
                     coloc_component = c(1L, 1L, 2L, NA_integer_), has_coloc = c(TRUE, TRUE, TRUE, FALSE))
b <- fig3_signal_counts(components, c("shared1", "shared2"))
pb <- fig3_panel_b(b)
expect(sum(b$signal_counts$nn) == 2L && sum(b$coloc_rate$n_cs) == 3L,
       "B preserves the distinct main and inset populations")
expect(inherits(pb, "wrapped_patch"), "B keeps one assembly tag around its inset")

# C and D refit the source model specification from a synthetic effect table.
set.seed(903)
effect <- tibble(posterior_mean.x = runif(120, -.3, .7),
                 posterior_mean.y = runif(120, -.3, .7),
                 log_maf = rnorm(120), max_af = rep(c(.001, .01, .05), 40))
models <- set_names(rep(list(effect), 3), c("eQTL-sQTL", "eQTL-pQTL", "sQTL-pQTL"))
c <- fig3_interaction_data(models)
pc <- fig3_panel_c(c)
expect(nrow(pc$data) == 3L && all(str_detect(pc$data$term, ":")),
       "C renders only the three newly fitted interaction terms")
reference <- broom::tidy(lm(log1p(posterior_mean.y) ~ log1p(posterior_mean.x) * log_maf,
                           data = effect), conf.int = TRUE) %>% filter(str_detect(term, ":"))
expect(abs(pc$data$estimate[pc$data$type == "eQTL-pQTL"] - reference$estimate) < 1e-12,
       "C preserves the signed eQTL-pQTL model")
d <- fig3_effect_data(models)
pd <- fig3_panel_d(d)
expect(identical(pd$data$maf_bin[1:3], c("Rare", "Low-Freq", "Common")),
       "D uses the source MAF bin boundaries")
expect(isTRUE(all.equal(pd$data$posterior_mean.x, log1p(effect$posterior_mean.x))),
       "D renders freshly transformed posterior effects")
expect(isTRUE(all.equal(pd$layers[[6]]$data$slope, d$slopes$slope)),
       "D labels the slopes refitted in each MAF bin")

# E filters fitted annotations, retaining both directions and source CI values.
enrichment <- tibble(feature = c("conserved", "PLS", "dELS", "intron_variant"),
                     fold = c(.6, 2, 1.3, .8), lo = c(.4, 1.2, 1.1, .5),
                     hi = c(.8, 3, 1.5, 1.2), p = c(.001, .002, .003, .4),
                     padj = c(.01, .02, .03, .8), n = c(10L, 8L, 4L, 30L))
pe <- fig3_panel_e(enrichment)
expect(identical(pe$data$feature, c("conserved", "PLS")),
       "E applies BH significance and the five-variant floor")
expect(identical(pe$data$dir, c("common", "rare")) && pe$coordinates$clip == "off",
       "E retains MAF direction and visible direction labels")
expect(isTRUE(all.equal(ggplot_build(pe)$data[[3]]$xmin, log10(pe$data$lo))),
       "E renders the fitted confidence limits on the log10 axis")
expect(inherits(fig3_panel_e(enrichment %>% mutate(padj = 1)), "ggplot"),
       "E renders an empty clean plot when no feature passes")

# Render the actual assembly to check all five panels and the B inset together.
figure <- fig3_assemble(list(A = pa, B = pb, C = pc, D = pd, E = pe))
grob <- patchwork::patchworkGrob(figure)
expect(inherits(grob, "gtable"), "the source-data assembly renders as five tagged panels")
tags <- character()
find_tags <- function(x) {
  if (inherits(x, "text")) {
    label <- paste(x$label, collapse = "")
    if (nchar(label) == 1L && label %in% LETTERS) tags <<- c(tags, label)
  } else if (inherits(x, "gtable")) {
    walk(x$grobs, find_tags)
  } else if (!is.null(x$children)) {
    walk(x$children, find_tags)
  }
  invisible(NULL)
}
find_tags(grob)
expect(identical(sort(tags), LETTERS[1:5]), "the rendered assembly has exactly tags A through E")
message("All inline Figure 3 checks passed.")
