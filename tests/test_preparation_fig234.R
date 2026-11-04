# Run from the repository root. Synthetic inputs contain no participant data.
suppressPackageStartupMessages(library(tidyverse))
source("preparation/Fig2/methods.R")
source("preparation/Fig2/current_methods.R")
source("preparation/Fig3/methods.R")
source("preparation/Fig3/symmetric_methods.R")
source("preparation/Fig3/current_methods.R")
source("preparation/Fig4/methods.R")
source("preparation/Fig4/current_methods.R")
source("preparation/Fig3/plot_bridge.R")
expect <- function(ok, label) {
  if (!isTRUE(ok)) stop("FAIL: ", label, call. = FALSE)
  message("PASS: ", label)
}
expect_error <- function(expr, label) {
  failed <- inherits(tryCatch(force(expr), error = identity), "error")
  expect(failed, label)
}

# Frequency thresholds are left-closed, and the unit is a variant–trait row.
freq <- tibble(variant = c("v", "v", "w", "x", "y", "z"), molecular_trait_id = letters[1:6],
               pip = c(0.91, 0.99, 0.91, 0.91, 0.9, 0.99),
               gvs_max_af = c(0.009, 0.01, 0.05, 0.1, 0.2, 0))
counts <- fig2_frequency_counts(list(Expression = freq))$counts
expect(identical(as.integer(counts$n), c(1L, 1L, 1L, 1L)), "Fig2A preserves pair rows and bin boundaries")
expect(abs(sum(counts$percent) - 100) < 1e-12, "Fig2A percentages sum to 100")
outlier <- tibble(focal_outliers = 0, focal_non_outliers = 8, bg_outliers = 3, bg_non_outliers = 9)
stat <- fig2_outlier_statistics(outlier)
expect(abs(stat$odds_ratio - (0.5 * 9.5) / (8.5 * 3.5)) < 1e-12, "Fig2E uses half-cell correction for OR")
expect(abs(stat$p_value - fisher.test(matrix(c(0, 8, 3, 9), 2, byrow = TRUE))$p.value) < 1e-12,
       "Fig2E Fisher test uses original counts")
cum <- suppressWarnings(make_cumulative_afc_curve_both(tibble(gene = c("g", "g", "g", "h"), log2_aFC = c(-1, -2, 1, 3))))
expect(all.equal(cum$cumulative_df %>% filter(gene == "g", direction == "down") %>% pull(expected_het_expr),
                 c(0.625, 0.5625)) == TRUE, "Fig2D orders and combines signed effects")
expect(!any(fig2_cumulative_display(cum)$plot_df$expected_het_expr >= 4), "Fig2D applies cap before display summaries")

# A common ALT allele in one ancestry must not become a rare MAF after max().
af <- as_tibble(setNames(rep(list(c(0.99, 0.9, NA_real_)), 7), AF_POP_COLS))
af$gvs_afr_af <- c(0.2, 0.1, 0.1)
expect(isTRUE(all.equal(.max_folded_maf(af), c(0.2, 0.1, NA_real_))), "Fig3 folds each AF before maximum and retains missingness")
features <- build_feature_matrix(tibble(max_af = c(0.001, 0.2), phylop = c(2.1, 2), distTSS = c(-1999, 2000),
                                        gene_type = "protein_coding", consequence = "intron_variant"))
expect(identical(features$conserved, c(TRUE, FALSE)) && identical(features$tss_proximal, c(TRUE, FALSE)),
       "Fig3 annotation cutoffs preserve the source phylop and signed-distance rules")
fdr <- add_fdr_threshold_label(tibble(p = c(0.95, 0.95, 0.1, NA)), "p", quiet = TRUE)
expect(identical(fdr$Coloc, c(TRUE, TRUE, FALSE, FALSE)), "Fig3 Bayesian FDR includes boundary ties and rejects missing probabilities")
rate <- coloc_cs_rate_fit(tibble(type = c("a", "a", "b", "b"), Coloc = c(TRUE, FALSE, TRUE, TRUE)))
expect(all(rate$conf.low >= 0 & rate$conf.high <= 1), "S22 Wilson limits remain in probability range")
expect(abs(coloc_cs_overall_rate(rate) - 0.75) < 1e-12, "S22 pooled rate uses counts")
u <- data.table::data.table(variant = c("v", "v", "w"), pip = c(0.25, 1, 0.8),
                           log_maf = c(-2, -2, -1), layer = c("eQTL", "pQTL", "pQTL"), flag = c(TRUE, FALSE, FALSE))
collapsed <- collapse_to_variants(u, "flag")
expect(abs(collapsed[variant == "v"]$pip - 0.5) < 1e-12 && nrow(collapsed) == 2L,
       "Fig3E collapses duplicate layers with geometric PIP weight")
set.seed(730)
v <- tibble(type = "eQTL-pQTL", Coloc = rep(c(TRUE, FALSE), 100),
            af = runif(200, 0.001, 0.4))
for (col in AF_POP_COLS) v[[col]] <- v$af
fit <- coloc_maf_rate_fit(v)
reference <- glm(Coloc ~ I(-log10(af)), family = binomial(), data = v)
expect(abs(log(fit$odds_ratio) - coef(reference)[2]) < 1e-10, "S23 uses unstandardized negative log10 MAF")

# Figure 4 counts unique genes, excludes failed artifact rows, and preserves keys.
gene <- tibble(gene = c("g", "g", "h"), trait = "trait", trait_category = "cat", any_coloc = TRUE)
rates <- tibble(trait = "trait", layer = "union", pct_stringent_union = 10, pct_lenient_union = 20)
trait <- build_plot_data(gene, rates, rate_type = "lenient")
expect(trait$count == 2L && trait$pct_coloc_union == 20, "Fig4A counts unique genes and selects lenient union rate")
expect_error(build_plot_data(gene, bind_rows(rates, rates)), "Fig4A rejects duplicate union-rate rows")
display_traits <- fig4_current_trait_data(tibble(trait_category = c(rep("category", 12), "neuropsychiatric"), count = 20:8))
expect(nrow(display_traits) == 10L && !any(display_traits$trait_category == "neuropsychiatric"),
       "Fig4A retains current category, count, and top-ten display filters")
tw <- tibble(phenotype_id = c("g", "g", "h", "h"), AF_threshold = c(0, 0.01, 0, 0.01),
             r.squared = c(0.2, 0.1, 0.3, 0.2), H2 = 0.5, SE = 0.1, Pval = c(1e-7, 1e-7, 1, 1))
r2 <- fig4_r2_data(tw)
expect(nrow(r2) == 1L && abs(r2$difference - 0.2) < 1e-12, "Fig4B filters heritability P and contrasts AF thresholds")
membership <- tibble(trait = c("t", "t", "t", "t", "u"), gene = c("g", "g", "g", "g", "h"),
                     layer = c("eQTL", "eQTL", "pQTL", "union", "sQTL"))
upset <- fig4_upset_counts(membership)
expect(sum(upset$intersections$gene_trait_pairs) == 2L && upset$multilayer_n == 1L,
       "S24 removes duplicate rows and excludes union layer")
qq <- calculate_qq_data(tibble(coloc_trait = "MCH", group = "test", common_pvalue = 0, p_rare_cond = 0.1), "MCH")
expect(all(is.finite(qq$observed)) && qq$observed_display[qq$test == "common"] == 40,
       "Fig4C floors zero P and applies display cap")
signals <- tibble(study = "study", trait = "monocyte", layer = "eQTL", gwas_cs = "gwas", qtl_sig = "ENSG00000101017:cs1",
                  gene = "ENSG00000101017", RCP = 0.9, CLPP = 0.2, RCP_pass_FDR = TRUE)
variants <- tibble(layer = "eQTL", qtl_cs = c("cs1", "unused"), position = c(46118278, 46118300),
                  variant = c("v1", "v2"), cs_id = c("cs1", "unused"), pip = c(0.99, 0.95), rsid = "rs", gene_id = "ENSG00000101017")
locus <- fig4_cd40_data(signals, variants)
expect(nrow(locus) == 1L && locus$variant == "v1" && locus$row_offset == 0,
       "S25 joins only participating credible sets and retains locus offsets")
reference <- tibble(assay = factor(c("a", "b")), estimate = c(1, 2))
rebuilt <- tibble(assay = c("a", "b"), estimate = c(1, 2))
expect(isTRUE(all.equal(reference, replace_prepared_frame(reference, rebuilt, "fixture"))),
       "Notebook bridge retains layout factor order and accepts recomputed values")
expect_error(replace_prepared_frame(reference, mutate(rebuilt, estimate = estimate + 1), "fixture"),
             "Notebook bridge stops when scientific values differ")
message("All Figure 2–4 preparation checks passed.")
