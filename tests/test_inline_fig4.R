# Synthetic checks for visible Figure 4 calculations and explicit plot layers.
suppressPackageStartupMessages({library(tidyverse); library(patchwork)})
source("tests/helpers.R")
source("preparation/Fig4/methods.R")
source("preparation/Fig4/current_methods.R")

plot_script <- "notebooks/current/plots/Fig4.R"
ok("Figure 4 has explicit plot construction", file.exists(plot_script))
source(plot_script)

# Changes to the heritability denominator or threshold contrast must fail here.
r2_source <- tibble(phenotype_id = rep(c("pass", "boundary"), each = 2),
                    AF_threshold = rep(c(0, 0.01), 2),
                    r.squared = c(0.3, 0.2, 0.9, 0.1), H2 = 0.5,
                    SE = 0.1, Pval = rep(c(1e-6, 1e-5), each = 2))
r2 <- fig4_r2_data(r2_source)
ok("R2 uses H2 and strict P threshold", nrow(r2) == 1L && abs(r2$difference - 0.2) < 1e-12)
r2_plot <- fig4_plot_r2(r2)
r2_bins <- suppressWarnings(ggplot_build(r2_plot))$data[[1]]
ok("R2 histogram counts the retained gene", sum(r2_bins$count) == 1L)
ok("R2 histogram preserves percent limits", identical(r2_plot$scales$get_scales("x")$limits, c(-20, 40)))

# Keep independent ranks within each group and test. Cap only display values.
matched <- tibble(coloc_trait = "MCH", group = rep(c("Colocalized genes", "Matched non-colocalized genes"), each = 2),
                  common_pvalue = c(0, 0.1, 0.01, 0.5), p_rare_cond = c(1e-20, 0.2, 0.05, 1))
qq <- calculate_qq_data(matched, "MCH") %>%
  mutate(group = if_else(str_detect(group, "Matched"), "Null", group))
ok("QQ retains uncapped evidence and caps display", max(qq$observed) > 40 && max(qq$observed_display[qq$test == "common"]) == 40 && max(qq$observed_display[qq$test == "rare"]) == 12.5)
qq_plot <- fig4_plot_qq(qq)
qq_points <- ggplot_build(qq_plot)$data[[2]]
ok("QQ layers use calculated ranks and values", isTRUE(all.equal(qq_points$x, qq$expected)) && isTRUE(all.equal(qq_points$y, qq$observed_display)))
ok("QQ retains both comparison colors", setequal(qq_points$colour, c("#275D8E", "grey50")))

enrichment <- tibble(type = factor(rep(c("Mendelian Gene\n(Any)", "Druggable gene"), each = 2)),
                     term = factor(rep(c("Rare", "Common"), 2), c("Rare", "Common")),
                     OR = c(1.2, 0.9, 1.3, 1.1), conf.low = c(1.05, 0.8, 1.15, 1.01),
                     conf.high = c(1.4, 1.05, 1.5, 1.2))
enrichment_plot <- fig4_plot_enrichment(enrichment)
intervals <- ggplot_build(enrichment_plot)$data[[2]]
ok("Enrichment plots supplied odds ratios and confidence limits", isTRUE(all.equal(intervals$x, enrichment$OR)) && isTRUE(all.equal(intervals$xmin, enrichment$conf.low)) && isTRUE(all.equal(intervals$xmax, enrichment$conf.high)))

traits <- tibble(trait = c("body mass index", "mchc"), trait_category = c("anthropometric", "hematological"),
                 count = c(20L, 10L), pct_coloc_union = c(25, 50))
trait_plot <- fig4_plot_traits(traits)
labels <- Filter(function(layer) inherits(layer$geom, "GeomText"), trait_plot$layers)[[1]]$data
ok("Trait labels retain current names and count-column position", setequal(labels$label, c("BMI", "MCHC")) && all(labels$x == 2000))
ok("Missing A leaves an assembly for the observed panels", inherits(fig4_assemble(NULL, r2_plot, qq_plot, enrichment_plot), "patchwork"))

# Duplicate membership records and the internal union layer cannot add pairs.
membership <- tibble(trait = c("t", "t", "t", "t", "t", "u"), gene = c("g", "g", "g", "g", "h", "g"),
                     layer = c("eQTL", "eQTL", "pQTL", "union", "sQTL", "sQTL"))
upset <- fig4_upset_counts(membership)
ok("S24 counts exclusive gene-trait pairs", identical(upset$intersections$gene_trait_pairs, c(2L, 1L)) && upset$multilayer_n == 1L && upset$all_three_n == 0L)

# One QTL variant can participate in two distinct GWAS loci. Each locus must
# retain its own offset, while unrelated credible sets and genes are excluded.
signals <- tibble(study = c("a", "b", "a"), trait = "monocyte", layer = "eQTL", gwas_cs = c("locus1", "locus2", "other"),
                  qtl_sig = "ENSG00000101017:cs1", gene = c("ENSG00000101017", "ENSG00000101017", "other"),
                  RCP = 0.9, CLPP = 0.2, RCP_pass_FDR = TRUE)
variants <- tibble(layer = "eQTL", qtl_cs = c("cs1", "unused"), position = c(46118278, 46118300),
                   variant = c("v1", "v2"), cs_id = c("cs1", "unused"), pip = 0.99, rsid = "rs", gene_id = "ENSG00000101017")
cd40 <- fig4_cd40_data(signals, variants)
ok("S25 preserves locus expansion and symmetric offsets", nrow(cd40) == 2L && all(cd40$variant == "v1") && isTRUE(all.equal(sort(cd40$row_offset), c(-0.065, 0.065))))
ok("S25 preserves genomic coordinates", all(cd40$position_mb == 46.118278) && all(cd40$trait_y == 1 + cd40$row_offset))

scripts <- c("notebooks/current/Fig4.R", "notebooks/current/supplements/S24.R", "notebooks/current/supplements/S25.R")
for (script in scripts) {
  text <- paste(readLines(script, warn = FALSE), collapse = "\n")
  ok(paste(script, "shows preparation methods"), grepl("# @include preparation/Fig4/current_methods.R", text, fixed = TRUE))
  ok(paste(script, "does not read saved results or call source"), !grepl("\\b(readRDS|source)\\s*\\(|objects/|current_panel_[A-D]", text, perl = TRUE))
}
ok("Plots have no titles or subtitles", all(map_lgl(list(r2_plot, qq_plot, enrichment_plot, trait_plot), ~is.null(.x$labels$title) && is.null(.x$labels$subtitle))))
message("Visible Figure 4 and S24/S25 checks passed.")
