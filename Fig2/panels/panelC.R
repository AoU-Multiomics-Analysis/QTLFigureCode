# Panel C — |log2 aFC| effect size by variant consequence (expression)
source("fig2/00_setup.R")
d <- readRDS(file.path(DATA_DIR, "C.rds"))

pC <- d$afc_plot_df %>%
  filter(!str_detect(color_group, "NA")) %>%
  ggplot(aes(abs_log2_aFC, name, color = color_group, fill = color_group)) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    stat_halfeye(fill_type = "segments", point_interval = "median_qi",
                 show.legend = FALSE, height = 1, point_size = 6,
                 interval_size = 0.4, slab_linewidth = 5.5, slab_alpha = 0.7) +
    scale_fill_manual(values = d$afc_palette) +
    scale_color_manual(values = d$afc_palette) +
    geom_text(data = d$log2_afc_regression, aes(x, name, label = label),
              inherit.aes = FALSE, hjust = 0, size = ANN_SM) +
    labs(x = "|Log2 aFC|", y = NULL) +
    coord_cartesian(clip = "off") +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          plot.margin = margin(5.5, 45, 5.5, 5.5))

saveRDS(pC, file.path(OBJ_DIR, "C.rds"))
message("Panel C built.")
