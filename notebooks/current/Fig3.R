# %% Read the current cross-QTL panels
#| The saved objects contain the current fitted results and plot data. This step
#| rebuilds their display. It does not rerun FastENLOC or the source coefficient tests.
#| The full upstream pipeline remains in Fig3/00_setup.R and Fig3/panels/.
suppressPackageStartupMessages({
  library(ggh4x)
  library(ggrepel)
})
panels <- set_names(LETTERS[1:5]) |>
  map(\(letter) readRDS(source_file(paste0("Fig3/objects/", letter, ".rds"))))

# Set text for the notebook canvas. The archived image retains its original headings.
notebook_theme <- theme(
  text = element_text(family = "sans", size = 8),
  axis.text = element_text(size = 7), axis.title = element_text(size = 8),
  strip.text = element_text(size = 8),
  plot.title = element_blank(), plot.subtitle = element_blank(),
  plot.margin = margin(3, 3, 3, 3, "mm"))
fit_panel <- function(plot) {
  before <- plot$data
  if (inherits(plot, "patchwork")) {
    plot$patches$plots <- map(plot$patches$plots, fit_panel)
    plot$patches$annotation$title <- NULL
    plot$patches$annotation$subtitle <- NULL
    plot$patches$annotation$theme <- theme()
  }
  for (i in seq_along(plot$layers)) {
    if (inherits(plot$layers[[i]]$geom, "GeomText")) {
      plot$layers[[i]]$aes_params$size <- min(c(plot$layers[[i]]$aes_params$size, 2.4))
    }
  }
  plot <- plot + notebook_theme + labs(title = NULL, subtitle = NULL)
  stopifnot(identical(before, plot$data))
  plot
}
panels <- map(panels, fit_panel)
panels$A$layers[[3]]$aes_params$size <- 2.1
for (i in seq_along(panels$D$layers)) {
  if (inherits(panels$D$layers[[i]]$geom, "GeomPoint")) panels$D$layers[[i]]$aes_params$size <- 0.9
}
# Keep the common/rare direction labels visible above panel E.
panels$E$coordinates$clip <- "off"

# %% Gene-level sharing
#| The denominator is 755 genes fine-mapped in all three assays. NumColocs
#| counts colocalized assay pairs. The percentages use this shared-gene set.
stopifnot(sum(panels$A$data$n) == 755L)
save_plot(panels$A, "Fig3/panel_A", 130, 90)
print(panels$A)

# %% Colocalized signal components
#| The saved patchwork includes the signal-count distribution and the per-assay
#| credible-set inset. The inset is distinct from the assay-pair rates in S22.
components <- readRDS(source_file("Fig3/data/coloc_cs_components.rds"))
signal_counts <- components |>
  filter(!is.na(coloc_component)) |>
  distinct(Gene, coloc_component) |>
  count(Gene) |>
  count(n, name = "n_genes") |>
  mutate(proportion = n_genes / sum(n_genes))
stopifnot(abs(sum(signal_counts$proportion) - 1) < 1e-12)
genes <- readRDS(source_file("Fig3/data/susie_finemapping.rds"))$genes_all_three_assays
rates <- components |> filter(Gene %in% genes) |> group_by(type) |>
  summarise(n_cs = n_distinct(cs_id), n_coloc = n_distinct(cs_id[has_coloc]),
            coloc_rate = n_coloc / n_cs, .groups = "drop")
# Compare both rebuilt summaries with the exact saved panel before restyling it.
saved_b <- attr(panels$B, "grobs")$full
stopifnot(isTRUE(all.equal(rates, saved_b$data, check.attributes = FALSE)),
          isTRUE(all.equal(signal_counts |> select(n, nn = n_genes),
                           saved_b$patches$plots[[1]]$data, check.attributes = FALSE)))
main <- ggplot(signal_counts, aes(factor(n), proportion, fill = factor(n))) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = c("#D8B7B7", "#C28F8F", "#B17373", "#A25B5B",
                               "#914747", "#803535", "#702626"), guide = "none") +
  scale_y_continuous(labels = scales::label_percent(), expand = expansion(mult = c(0, .05))) +
  labs(x = "Number of colocalized signals", y = "Percent of genes") + notebook_theme
inset <- ggplot(rates, aes(coloc_rate, reorder(type, coloc_rate), fill = type)) +
  geom_col(show.legend = FALSE) +
  scale_fill_manual(values = c(eQTL = "#007C8D", pQTL = "#D44258", sQTL = "#EFB144")) +
  scale_x_continuous(labels = scales::label_percent(), expand = expansion(mult = c(0, .05))) +
  labs(x = "% CS coloc", y = NULL) + notebook_theme +
  theme(axis.text = element_text(size = 6), axis.title = element_text(size = 7),
        plot.margin = margin(0, 0, 0, 0))
panels$B <- wrap_elements(full = main + inset_element(
  inset, left = .48, bottom = .60, right = .98, top = .98, ignore_tag = TRUE))
save_plot(panels$B, "Fig3/panel_B", 150, 100)
print(panels$B)

# %% Effect-size concordance and allele frequency
#| Panel C uses the saved interaction coefficients and 95% confidence limits.
#| Panel D uses the saved transformed expression and protein effects by MAF.
#| These are the recorded models; no new model specification is substituted.
stopifnot(nrow(panels$C$data) == 3L, nrow(panels$D$data) == 1722L)
save_plot(panels$C, "Fig3/panel_C", 140, 100)
print(panels$C)
save_plot(panels$D, "Fig3/panel_D", 180, 100)
print(panels$D)

# %% Variant annotations
#| Panel E contains the current annotation results. It is no longer the blank
#| placeholder described by the original assembly-script comment.
stopifnot(nrow(panels$E$data) > 0L)
save_plot(panels$E, "Fig3/panel_E", 150, 100)
print(panels$E)

# %% Assemble the current five panels
#| Keep the current panel order. The exact manuscript image below also contains
#| the final vector lettering and the archived export resolution.
coloc_counts <- panels$A + panels$B + plot_layout(widths = c(1, 1))
top_row <- (coloc_counts | panels$C) + plot_layout(widths = c(2.5, 1))
bottom_row <- (panels$D | panels$E) + plot_layout(widths = c(1.3, 1))
fig3 <- (top_row / bottom_row) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 11, face = "bold"))
save_plot(fig3, "Fig3/Figure3_rebuilt", 184, 130)
print(fig3)

# %% Exact submission reference
#| This archived PNG is the Figure 3 image embedded in the November 2 manuscript.
show_reference(3)
