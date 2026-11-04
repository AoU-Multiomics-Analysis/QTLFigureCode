# %% Read source credible sets for all populations
#| Aggregate variant rows to one credible set and one lead PIP per set. Include the current EAS eQTL input.
# @include notebooks/current/analysis/Fig1_helpers.R
fig1_groups <- c("COMB", "AFR", "AMR", "EUR", "SAS", "EAS", "MID")
# @include notebooks/current/analysis/Fig1_finemapping.R
labels <- c(COMB = "Combined", AFR = "AFR", AMR = "AMR", EUR = "EUR", SAS = "SAS", EAS = "EAS", MID = "MID")
lead <- sets %>% select(cs_id, lead_pip, lead_variant, n_tied_lead_variants, assay, population)

# %% Summarize lead PIPs and estimate reflected densities
#| Lead PIP is the maximum PIP within each credible set. Each set has equal weight.
#| Use bandwidth 0.04, reflect the density at 0 and 1, and integrate the displayed curve as a check.
summary <- lead |>
  group_by(assay, population) |>
  summarize(n_credible_sets = n(), mean_lead_pip = mean(lead_pip),
            median_lead_pip = median(lead_pip), min_lead_pip = min(lead_pip),
            max_lead_pip = max(lead_pip),
            percent_pip_above_09 = 100 * mean(lead_pip > .9),
            sets_with_tied_leads = sum(n_tied_lead_variants > 1), .groups = 'drop')
order_check <- summary |> arrange(assay, population) |> group_by(assay) |>
  summarize(decreasing_mean_PIP = all(diff(mean_lead_pip) <= 0),
            decreasing_resolved_percent = all(diff(percent_pip_above_09) <= 0), .groups = 'drop')
stopifnot(all(order_check$decreasing_mean_PIP), all(order_check$decreasing_resolved_percent))

colours <- c(COMB = '#822E4A', AFR = '#6F7F20', AMR = '#39743F', EAS = '#9B5900',
             EUR = '#06664C', SAS = '#917502', MID = '#9A3B2D')
types <- c(COMB = 'solid', AFR = 'solid', AMR = '22', EAS = '13', EUR = 'solid', SAS = '44', MID = 'solid')
bw <- .04
estimator <- ggplot(lead, aes(lead_pip, group = population)) +
  geom_density(bw = bw, bounds = c(0, 1), n = 2048) + facet_grid(~assay) +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = 0))
density_data <- ggplot_build(estimator)$data[[1]]
curves <- density_data |>
  transmute(assay = factor(assays[as.integer(PANEL)], levels = assays),
            population = factor(groups[group], levels = groups), pip = x,
            density = density, n_sets = n)
area <- curves |>
  group_by(assay, population) |> arrange(pip, .by_group = TRUE) |>
  summarize(area = sum(diff(pip) * (head(density, -1) + tail(density, -1)) / 2),
            n_sets = first(n_sets), .groups = 'drop') |>
  left_join(summary |> select(assay, population, n_credible_sets), by = c('assay', 'population'))
stopifnot(nrow(area) == 21L, all(abs(area$area - 1) < .01),
          all(area$n_sets == area$n_credible_sets), all(curves$density >= 0))
# %% Draw the density, mean PIP, and high-PIP percentage panels
#| The percentage uses the strict threshold PIP > 0.9. Groups can contain different detected signals.
base_theme <- theme_classic(base_size = 8, base_family = 'Arial') +
  theme(axis.text = element_text(size = 7.5, colour = '#263B45'),
        axis.title = element_text(size = 8.5),
        axis.title.x = element_text(margin = margin(t = 4)),
        axis.title.y = element_text(margin = margin(r = 4)),
        axis.line = element_line(linewidth = .25, colour = '#53616A'),
        axis.ticks = element_line(linewidth = .25, colour = '#53616A'),
        axis.ticks.length = unit(.8, 'mm'),
        strip.background = element_blank(),
        strip.text = element_text(size = 9, face = 'bold', margin = margin(b = 4)),
        panel.spacing.x = unit(6, 'mm'), plot.margin = margin(4, 2, 2, 2, 'mm'))
a <- ggplot(curves, aes(pip, density, colour = population, linetype = population)) +
  geom_vline(xintercept = .9, colour = '#ACB5BA', linewidth = .25, linetype = '22') +
  geom_line(linewidth = .6, key_glyph = 'path') + facet_grid(~assay) +
  scale_colour_manual(values = colours, breaks = groups, labels = labels, name = NULL) +
  scale_linetype_manual(values = types, breaks = groups, labels = labels, name = NULL) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, .25, .5, .75, 1),
                     labels = c('0', '0.25', '0.5', '0.75', '1'), expand = expansion(mult = 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, .07))) +
  labs(x = 'Lead-variant PIP', y = 'Density') + base_theme +
  theme(legend.position = 'bottom', legend.text = element_text(size = 7.5),
        legend.key.width = unit(6, 'mm'), legend.key.height = unit(3, 'mm'),
        legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, 0, 0)) +
  guides(colour = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1))
bdata <- summary |> mutate(population = factor(population, levels = rev(groups)),
                           label = sprintf('%.3f', mean_lead_pip))
b <- ggplot(bdata, aes(mean_lead_pip, population, colour = population)) +
  geom_point(aes(shape = population), size = 2.2, show.legend = FALSE) +
  geom_text(aes(x = mean_lead_pip + .035, label = label), hjust = 0,
            size = 7.5 / .pt, family = 'Arial', colour = '#263B45', show.legend = FALSE) +
  facet_grid(~assay) + scale_colour_manual(values = colours, guide = 'none') +
  scale_shape_manual(values = c(COMB = 18, AFR = 16, AMR = 16, EAS = 16, EUR = 16, SAS = 16, MID = 16), guide = 'none') +
  scale_y_discrete(labels = labels, expand = expansion(add = .6)) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, .25, .5, .75, 1),
                     labels = c('0', '0.25', '0.5', '0.75', '1'), expand = expansion(mult = 0)) +
  labs(x = 'Mean lead-variant PIP', y = NULL) + base_theme +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())
cdata <- summary |>
  mutate(population = factor(population, levels = rev(groups)),
         label = sprintf('%.1f%%', percent_pip_above_09))
fills <- c(COMB = '#CF597E', AFR = '#C8D892', AMR = '#7FC385', EAS = '#ECA676',
           EUR = '#30AC86', SAS = '#EBCE8D', MID = '#E47D72')
c <- ggplot(cdata, aes(percent_pip_above_09, population,
                       colour = population, fill = population)) +
  geom_col(width = .62, linewidth = .25, show.legend = FALSE) +
  geom_text(aes(x = percent_pip_above_09 + 2, label = label), hjust = 0,
            size = 7.5 / .pt, family = 'Arial', colour = '#263B45', show.legend = FALSE) +
  facet_grid(~assay) + scale_colour_manual(values = colours, guide = 'none') +
  scale_fill_manual(values = fills, guide = 'none') +
  scale_y_discrete(labels = labels, expand = expansion(add = .6)) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100),
                     labels = function(x) paste0(x, '%'), expand = expansion(mult = 0)) +
  labs(x = 'Credible sets with lead-variant PIP > 0.9', y = NULL) + base_theme +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())
tag_theme <- theme(plot.tag = element_text(size = 11, face = 'bold', family = 'Arial'),
                   plot.tag.position = c(0, 1))
p <- ((a + tag_theme) / (b + tag_theme) / (c + tag_theme)) +
  plot_layout(heights = c(1.2, 1, 1)) + plot_annotation(tag_levels = 'A')
stopifnot(nrow(ggplot_build(b)$data[[1]]) == 21L,
          max(abs(sort(ggplot_build(b)$data[[1]]$x) - sort(summary$mean_lead_pip))) < 1e-12,
          nrow(ggplot_build(c)$data[[1]]) == 21L,
          max(abs(sort(ggplot_build(c)$data[[1]]$x) - sort(summary$percent_pip_above_09))) < 1e-12)

figure <- p
print(figure)

# %% Save aggregate outputs
save_plot(figure, "supplements/S13", 184, 210)
write_tsv(summary, output_file("supplements/S13_summary.tsv"))
