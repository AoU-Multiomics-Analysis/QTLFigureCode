# Internal helper for 05_supplements.R. Input tables come from 04_finemapping.R.
# EAS update of the existing supplemental analysis; no model is refitted.
suppressPackageStartupMessages({library(tidyverse); library(patchwork); library(svglite); library(jsonlite)})
root <- file.path(config$output_root, 'paper')
base <- file.path(root, 'output/eas_fine_mapping_2026-11-02')
assays <- c('eQTL', 'sQTL', 'pQTL')
groups <- c('COMB', 'AFR', 'AMR', 'EUR', 'SAS', 'EAS', 'MID')
labels <- c(COMB = 'Combined', AFR = 'AFR', AMR = 'AMR', EAS = 'EAS', EUR = 'EUR', SAS = 'SAS', MID = 'MID')
cs <- read_tsv(file.path(base, 'data/credible_sets.tsv.gz'), show_col_types = FALSE) |>
  mutate(assay = factor(assay, levels = assays), population = factor(population, levels = groups))
out <- file.path(base, 'credible_set_sizes')
dir.create(out, showWarnings = FALSE)
summary <- cs |>
  group_by(assay, population) |>
  summarize(n_cs = n(), mean_size = mean(cs_size),
            min_size = min(cs_size), q1 = quantile(cs_size, .25),
            median_size = median(cs_size), q3 = quantile(cs_size, .75),
            max_size = max(cs_size),
            log_q1 = quantile(log10(cs_size), .25),
            log_median = median(log10(cs_size)),
            log_q3 = quantile(log10(cs_size), .75), .groups = "drop") |>
  mutate(log_fence_low = log_q1-1.5*(log_q3-log_q1),
         log_fence_high = log_q3+1.5*(log_q3-log_q1))
whiskers <- cs |>
  left_join(summary |> select(assay, population, log_fence_low, log_fence_high),
            by = c("assay", "population")) |>
  mutate(outlier = log10(cs_size) < log_fence_low | log10(cs_size) > log_fence_high) |>
  group_by(assay, population) |>
  summarize(whisker_low = min(cs_size[!outlier]),
            whisker_high = max(cs_size[!outlier]), n_outliers = sum(outlier), .groups = "drop")
summary <- summary |> left_join(whiskers, by = c("assay", "population"))
order_check <- summary |> arrange(assay, population) |> group_by(assay) |>
  summarize(increasing_mean_size = all(diff(mean_size) >= 0), .groups = 'drop')
stopifnot(all(order_check$increasing_mean_size))
fills <- c(COMB = '#CF597E', AFR = '#C8D892', AMR = '#7FC385', EAS = '#ECA676', EUR = '#30AC86', SAS = '#EBCE8D', MID = '#E47D72')
lines <- c(COMB = '#822E4A', AFR = '#6F7F20', AMR = '#39743F', EAS = '#9B5900', EUR = '#06664C', SAS = '#917502', MID = '#9A3B2D')
p <- ggplot(cs, aes(population, cs_size, fill = population, colour = population)) +
  geom_boxplot(width = .58, linewidth = .35, outlier.size = .45,
               outlier.stroke = .1, outlier.alpha = .8, show.legend = FALSE) +
  facet_wrap(vars(assay), ncol = 1, axes = "all", axis.labels = "all") +
  scale_fill_manual(values = fills, guide = "none") +
  scale_colour_manual(values = lines, guide = "none") +
  scale_x_discrete(limits = groups,
                   labels = c(COMB = "Combined", AFR = "AFR", AMR = "AMR", EAS = "EAS", EUR = "EUR", SAS = "SAS", MID = "MID"),
                   expand = expansion(add = .65)) +
  scale_y_log10(limits = c(.8, 3300), breaks = c(1, 10, 100, 1000),
                labels = scales::label_comma(), expand = expansion(mult = 0)) +
  labs(x = "Population / analysis", y = "Credible set size (variants)") +
  theme_classic(base_size = 8, base_family = "Arial") +
  theme(axis.text = element_text(size = 7.5, colour = "#263B45"),
        axis.text.x = element_text(margin = margin(t = 3)),
        axis.title = element_text(size = 8),
        axis.title.x = element_text(margin = margin(t = 5)),
        axis.title.y = element_text(margin = margin(r = 5)),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.y = element_line(linewidth = .3, colour = "#263B45"),
        axis.ticks.y = element_line(linewidth = .25, colour = "#263B45"),
        axis.ticks.length = unit(.8, "mm"),
        strip.background = element_blank(),
        strip.text = element_text(size = 9, face = "bold", margin = margin(b = 4)),
        panel.spacing.y = unit(4, "mm"),
        plot.margin = margin(3, 3, 3, 3, "mm"))
built <- ggplot_build(p)
observed <- built$data[[1]] |>
  transmute(assay = factor(assays[as.integer(PANEL)], levels = assays),
            population = factor(groups[group], levels = groups),
            observed_q1 = lower, observed_median = middle, observed_q3 = upper,
            observed_low = ymin, observed_high = ymax,
            observed_n_outliers = map_int(outliers, length))
check <- summary |> left_join(observed, by = c("assay", "population"))
stopifnot(nrow(check) == 21L,
          all(abs(check$log_q1-check$observed_q1) < 1e-10),
          all(abs(check$log_median-check$observed_median) < 1e-10),
          all(abs(check$log_q3-check$observed_q3) < 1e-10),
          all(abs(log10(check$whisker_low)-check$observed_low) < 1e-10),
          all(abs(log10(check$whisker_high)-check$observed_high) < 1e-10),
          all(check$n_outliers == check$observed_n_outliers))


ref <- tibble(assay = factor('pQTL', levels = assays),
              population = factor('EUR', levels = groups), reference_size = 9.7, label = 'UKBB pQTL')
eur_position <- match('EUR', groups)
p <- p +
  geom_segment(data = ref, aes(x = !!(eur_position + .03), xend = !!(eur_position + .24), y = 10.8, yend = 18),
               inherit.aes = FALSE, colour = 'black', linewidth = .25) +
  geom_point(data = ref, aes(population, reference_size), inherit.aes = FALSE,
             shape = 23, size = 2.1, stroke = .35, fill = '#4B9FC3', colour = '#186A91') +
  geom_text(data = ref, aes(x = !!(eur_position + .28), y = 22, label = label), inherit.aes = FALSE,
            family = 'Arial', fontface = 'bold', size = 7.5 / .pt, hjust = 0, colour = 'black')
saved <- list(plot = p)
write_tsv(summary, file.path(out, 'population_credible_set_summary.tsv'))
write_tsv(cs, file.path(out, 'credible_sets.tsv.gz'))
write_tsv(ref, file.path(out, 'UKBB_reference_marker.tsv'))
cs <- cs |> mutate(gene_id = str_extract(cs_id, 'ENSG[0-9]+'))
stopifnot(!anyNA(cs$gene_id))
gene_sums <- cs |>
  group_by(assay, population, gene_id) |>
  summarize(n_cs = n(), sum_size = sum(cs_size), .groups = "drop")

# Estimate the difference between two ratio means using gene-level influence
# contributions. The AFR contribution is subtracted within the same gene.
# This includes uncertainty in AFR and covariance from shared genes. No
# bootstrap, permutation, or fine-mapping refit is performed.
contrasts <- map_dfr(assays, function(a) {
  message("Calculate analytic differences and SEs: ", a)
  g <- gene_sums |> filter(assay == a)
  genes <- sort(unique(g$gene_id))
  G <- length(genes)
  N <- matrix(0, nrow = G, ncol = length(groups), dimnames = list(genes, groups))
  S <- N
  at <- cbind(match(g$gene_id, genes), match(as.character(g$population), groups))
  N[at] <- g$n_cs
  S[at] <- g$sum_size
  n_total <- colSums(N)
  mu <- colSums(S)/n_total
  residual <- S-sweep(N, 2, mu, `*`)
  contribution <- sweep(residual, 2, n_total, `/`)
  stopifnot(max(abs(colSums(contribution))) < 1e-10)
  covariance <- (G/(G-1))*crossprod(contribution)
  afr <- match("AFR", groups)
  delta_contribution <- sweep(contribution, 1, contribution[, afr], `-`)
  variance_delta <- (G/(G-1))*colSums(delta_contribution^2)
  # Independently express the same variance as V(group)+V(AFR)-2*Cov.
  second_variance <- diag(covariance)+covariance[afr, afr]-2*covariance[, afr]
  stopifnot(max(abs(variance_delta-second_variance)) < 1e-10)
  tibble(assay = a, population = groups, n_cs = n_total,
         n_cs_AFR = n_total[[afr]], mean_size = mu,
         AFR_mean_size = mu[[afr]], difference = mu-mu[[afr]],
         standard_error = sqrt(pmax(variance_delta, 0)),
         group_mean_SE = sqrt(diag(covariance)),
         AFR_mean_SE = sqrt(covariance[afr, afr]),
         covariance_with_AFR = covariance[, afr],
         gene_clusters_in_assay = G, n_genes_group = colSums(N > 0),
         n_genes_AFR = sum(N[, afr] > 0),
         n_shared_genes_with_AFR = colSums((N > 0) & (N[, afr] > 0)))
}) |>
  mutate(assay = factor(assay, levels = assays),
         population = factor(population, levels = groups),
         lower_1SE = difference-standard_error,
         upper_1SE = difference+standard_error)
expected <- cs |>
  group_by(assay, population) |>
  summarize(expected_mean = mean(cs_size), expected_n = n(), .groups = "drop")
checks <- contrasts |> left_join(expected, by = c("assay", "population"))
stopifnot(all(abs(checks$mean_size-checks$expected_mean) < 1e-10),
          all(checks$n_cs == checks$expected_n),
          all(contrasts$standard_error >= 0),
          all(contrasts$difference[contrasts$population == "AFR"] == 0),
          all(contrasts$standard_error[contrasts$population == "AFR"] == 0))

fills <- c(COMB = '#CF597E', AFR = '#C8D892', AMR = '#7FC385', EAS = '#ECA676', EUR = '#30AC86', SAS = '#EBCE8D', MID = '#E47D72')
lines <- c(COMB = '#822E4A', AFR = '#6F7F20', AMR = '#39743F', EAS = '#9B5900', EUR = '#06664C', SAS = '#917502', MID = '#9A3B2D')
bdata <- contrasts |>
  filter(population != "AFR") |>
  mutate(population = factor(population, levels = rev(setdiff(groups, 'AFR'))))
plot_b <- ggplot(bdata, aes(y = population, x = difference, colour = population, fill = population)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "#70818A", linewidth = .3) +
  geom_errorbar(aes(xmin = lower_1SE, xmax = upper_1SE),
                 orientation = "y", width = .32, linewidth = .4, show.legend = FALSE) +
  geom_point(shape = 21, size = 1.6, stroke = .35, show.legend = FALSE) +
  facet_grid(~ assay) +
  scale_colour_manual(values = lines, guide = "none") +
  scale_fill_manual(values = fills, guide = "none") +
  scale_y_discrete(labels = c(COMB = "Combined", AMR = "AMR", EAS = "EAS", EUR = "EUR", SAS = "SAS", MID = "MID"),
                   expand = expansion(add = .55)) +
  scale_x_continuous(limits = c(-5, 20), breaks = seq(-5, 20, 5),
                     expand = expansion(mult = 0)) +
  labs(x = "Mean credible-set size difference vs AFR (variants)", y = NULL) +
  theme_classic(base_size = 8, base_family = "Arial") +
  theme(axis.text = element_text(size = 7.5, colour = "#263B45"),
        axis.text.x = element_text(margin = margin(t = 3)),
        axis.text.y = element_text(margin = margin(r = 3)),
        axis.title.x = element_text(size = 8, margin = margin(t = 5)),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.x = element_line(linewidth = .3, colour = "#263B45"),
        axis.ticks.x = element_line(linewidth = .25, colour = "#263B45"),
        axis.ticks.length = unit(.8, "mm"),
        strip.background = element_blank(),
        strip.text = element_text(size = 9, face = "bold", margin = margin(b = 5)),
        panel.spacing.x = unit(4, "mm"),
        plot.margin = margin(3, 3, 3, 3, "mm"))
built_b <- ggplot_build(plot_b)
stopifnot(nrow(built_b$data[[2]]) == 18L, nrow(built_b$data[[3]]) == 18L,
          all(bdata$lower_1SE > -5), all(bdata$upper_1SE < 20),
          !anyNA(built_b$data[[2]]$xmin), !anyNA(built_b$data[[2]]$xmax))
plot_a <- saved$plot
stopifnot(identical(ggplot_build(plot_a)$data, ggplot_build(saved$plot)$data))
render <- function(p, name, height) {
  ggsave(file.path(out, paste0(name, ".pdf")), p, width = 184, height = height,
         units = "mm", device = cairo_pdf, bg = "white")
  ggsave(file.path(out, paste0(name, ".svg")), p, width = 184, height = height,
         units = "mm", device = svglite::svglite, bg = "white")
}
render(plot_a, "panel_A_native", 124)
render(plot_b, "panel_B_native", 63)
write_tsv(contrasts, file.path(out, "AFR_mean_size_differences.tsv"))
write_tsv(gene_sums, file.path(out, "gene_cluster_summaries.tsv.gz"))

write_lines(c(
 'Credible-set sizes and differences from AFR across the combined analysis and six populations.',
 '(A) Boxes show the median and interquartile range of log10-transformed credible-set sizes.',
 'Whiskers extend to the most extreme value within 1.5 interquartile ranges. Points show outliers.',
 'All observed sizes are retained. The blue diamond at EUR shows the UKBB reference mean of 9.7 variants.',
 '(B) Points show each group mean minus the AFR mean, in variants. Bars show plus or minus one standard error.',
 'The dashed line marks the AFR reference of zero. Negative values indicate smaller mean sets than AFR.',
 'Standard errors use stable Ensembl genes as clusters and include covariance between groups.',
 'For gene i and group g, z_ig = (S_ig - mean_g*N_ig)/N_g.',
 'S_ig is the sum of sizes. N_ig is the gene/group set count. N_g is the total group set count.',
 'SE(group-AFR) = sqrt[G/(G-1) * sum_i (z_ig-z_iAFR)^2]. G is the union of genes within each assay.',
 'Combined is a separate fine-mapping analysis. EAS means East Asian.',
 'Groups can contain different detected signals. Independent component labels do not identify matched signals.',
 'No fine-mapping model was refitted. No resampling was performed.'), file.path(out, 'figure_caption.txt'))
saveRDS(list(A = plot_a, B = plot_b, differences = contrasts), file.path(out, 'figure_display.rds'))
jsonlite::write_json(list(total_credible_sets = nrow(cs), groups = groups, EAS_included = TRUE,
                UKBB_reference_retained = TRUE, legacy_set_sizes_unchanged = TRUE,
                gene_cluster_key = 'stable Ensembl gene ID without version',
                covariance_with_AFR_included = TRUE, new_models_run = FALSE,
                standard_error_formula = 'sqrt[G/(G-1) * sum_i(z_ig-z_iAFR)^2]'),
           file.path(out, 'build_audit.json'), auto_unbox = TRUE, pretty = TRUE)
