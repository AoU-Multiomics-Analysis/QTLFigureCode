# %% Read source credible sets for all populations
#| Aggregate variant rows to one credible set and one lead PIP per set. Include the current EAS eQTL input.
# @include notebooks/current/analysis/Fig1_helpers.R
fig1_groups <- c("COMB", "AFR", "AMR", "EUR", "SAS", "EAS", "MID")
# @include notebooks/current/analysis/Fig1_finemapping.R
labels <- c(COMB = "Combined", AFR = "AFR", AMR = "AMR", EUR = "EUR", SAS = "SAS", EAS = "EAS", MID = "MID")
cs <- sets %>% select(assay, population, cs_id, cs_size)

# %% Calculate size distributions and draw their boxplots
#| The boxplot uses log10 sizes, retains all outliers, and shows the supplied UKBB pQTL reference mean of 9.7 variants.
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

# %% Calculate gene-cluster uncertainty for differences from AFR
#| Each gene contributes its set count and size sum. Subtract the AFR influence within each gene to include covariance.
cs <- cs |> mutate(gene_id = str_extract(cs_id, 'ENSG[0-9]+'))

stopifnot(!anyNA(cs$gene_id))

gene_sums <- cs |>
  group_by(assay, population, gene_id) |>
  summarize(n_cs = n(), sum_size = sum(cs_size), .groups = "drop")

contrasts <- fig1_size_contrasts(cs, assays, groups)

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

plot_a <- p


# %% Assemble the size and contrast panels
figure <- plot_a / plot_b + plot_layout(heights = c(124, 63)) + plot_annotation(tag_levels = "A")
print(figure)

# %% Save aggregate outputs
save_plot(figure, "supplements/S12", 184, 195)
write_tsv(cs %>% count(assay, population, cs_size, name = "credible_sets"),
          output_file("supplements/S12_size_frequencies.tsv"))
write_tsv(contrasts, output_file("supplements/S12_AFR_contrasts.tsv"))
