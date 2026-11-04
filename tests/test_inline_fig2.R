# These fixtures test the plots made from source-result summaries.
suppressPackageStartupMessages(library(tidyverse))
source("tests/helpers.R")
source("preparation/Fig2/methods.R")
source("preparation/Fig2/current_methods.R")

plot_file <- "notebooks/current/plots/Fig2.R"
ok("Figure 2 can construct plots from fresh results", file.exists(plot_file))
plot_env <- new.env(parent = globalenv())
plot_env$readRDS <- function(...) stop("A plot must not read a saved object")
sys.source(plot_file, envir = plot_env)

# Dropping duplicate variant IDs would change the rare bin from 40% to 25%.
freq <- tibble(variant = c("v", "v", "w", "x", "y", "excluded"),
               molecular_trait_id = letters[1:6], pip = c(rep(.99, 5), .9),
               gvs_max_af = c(.001, .002, .01, .05, .10, .001))
counts <- fig2_frequency_counts(list(Expression = freq))$counts
a <- plot_env$fig2_plot_a(counts)
ok("A counts variant-trait rows in the four current bins",
   identical(as.numeric(a$data$percent), c(40, 20, 20, 20)))
ok("A retains the 1%, 5%, and 10% boundaries",
   identical(as.character(a$data$frequency_bin), c("<1%", "1–<5%", "5–<10%", "10–50%")))

# The central reference bin must not appear in E. Its empty position remains in F.
outliers <- tibble(focal_lower_bound = rep(c(-100, -10, 175), 2),
                  focal_outliers = c(0, 2, 3, 4, 2, 0), focal_non_outliers = 8,
                  bg_outliers = 3, bg_non_outliers = 9,
                  type = rep(c("Down", "Up"), each = 3))
outliers <- fig2_outlier_statistics(outliers)
median <- tibble(pct_bin = c("(-100,-10)", "(-10,10)", "(175,200)"),
                 lower = c(-100, -10, 175), gene_type = "Protein Coding",
                 median_genes = c(2, 100, 3), q25_genes = c(1, 90, 2), q75_genes = c(3, 110, 4))
e <- plot_env$fig2_plot_e(outliers, median)
ok("E omits the reference bin and retains corrected odds ratios",
   nrow(e$data) == 4L && !any(e$data$focal_lower_bound == -10) &&
     abs(e$data$odds_ratio[1] - (.5 * 9.5) / (8.5 * 3.5)) < 1e-12)
ok("E maps bin centers to the same category positions as F",
   isTRUE(all.equal(plot_env$fig2_bin_transform(median$lower)$transform(c(-100, -10, 175)), c(1, 2, 3))))

expr <- tibble(variant = rep("same_variant", 60), pip = .99,
               gvs_all_af = seq(.001, .49, length.out = 60), gvs_max_af = gvs_all_af,
               gvs_max_subpop = rep(c("afr", "oth"), 30),
               gvs_odds_ratio = 10, gvs_p_value = .0001)
for (pop in c("afr", "amr", "eas", "eur", "mid", "sas", "oth")) expr[[paste0("gvs_", pop, "_af")]] <- .01
expr$gvs_oth_af[seq(1, 60, 2)] <- 0
bins <- fig2_ancestry_bins(expr)
s16 <- plot_env$fig2_plot_s16(bins)
ok("S16 keeps the row-level BH universe when variant IDs repeat",
   sum(bins$stack_df$n_enriched) == 60 && sum(distinct(bins$stack_df, af_bin, total_variants)$total_variants) == 60)
ok("S16 checks all seven populations for the shared-variant curve",
   all(bins$line_df$prop_shared == .5) && all(bins$line_df$prop_shared_smooth == .5))
fill <- ggplot_build(s16)$plot$scales$get_scales("fill")
ok("S16 displays OTH as REM while retaining source ancestry keys",
   "OTH" %in% s16$data$pop && identical(fill$get_labels(), c("AFR", "REM")))
ok("S16 retains the current ancestry colors shared with Figure 1",
   identical(unname(fill$map(c("AFR", "OTH"))), c("#C8D892", "#009392")))
missing_af <- expr
missing_af$gvs_all_af[1] <- NA_real_
missing_plot <- plot_env$fig2_plot_s16(suppressWarnings(fig2_ancestry_bins(missing_af)))
plot_warnings <- character()
invisible(suppressMessages(withCallingHandlers(ggplot_build(missing_plot), warning = function(w) {
  plot_warnings <<- c(plot_warnings, conditionMessage(w))
  invokeRestart("muffleWarning")
})))
ok("S16 retains missing-AF summaries without plot-removal warnings",
   anyNA(missing_plot$data$x) && !length(plot_warnings))
message("All inline Figure 2 checks passed.")
