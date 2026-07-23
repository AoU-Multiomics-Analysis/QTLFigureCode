# Panel B — signal-level colocalization: percent of genes by the number of
# colocalized signals, with an inset of per-assay credible-set coloc rate.
# Inputs: enloc_cs, susie_eqtl, susie_sqtl, susie_pqtl
source("fig3/00_setup.R")

cs_level        <- coloc_cs_components()
genes_all_three_assays <- susie_finemapping()$genes_all_three_assays

signal_counts <- cs_level %>%
  filter(!is.na(coloc_component)) %>%       # the original used na.omit()
  distinct(Gene, coloc_component) %>%
  dplyr::count(Gene) %>%
  dplyr::count(n, name = "nn")

main <- signal_counts %>%
  ggplot(aes(x = factor(n), y = nn / sum(nn), fill = factor(n))) +
    geom_col(width = 0.75) +
    scale_fill_manual(values = SIGNAL_FILL, guide = "none") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, 0.55)) +
    labs(x = "Number of colocalized signals", y = "Percent of genes",
         title = "Signal-level Colocalization") +
    theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), axis.line.y = element_blank(),
          plot.title = element_text(hjust = 0.5, size = SZ_TITLE, face = "bold"),
          plot.title.position = "plot")

coloc_rate <- cs_level %>%
  filter(Gene %in% genes_all_three_assays) %>%
  group_by(type) %>%
  summarise(n_cs       = n_distinct(cs_id),
            n_coloc    = n_distinct(cs_id[has_coloc]),
            coloc_rate = n_coloc / n_cs,
            .groups = "drop")

inset <- coloc_rate %>%
  ggplot(aes(y = reorder(type, coloc_rate), x = coloc_rate, fill = type)) +
    geom_col(show.legend = FALSE) +
    paletteer::scale_fill_paletteer_d("ltc::minou") +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                       expand = c(0, 0)) +
    labs(x = "% CS coloc", y = NULL) +
    theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
          axis.line.y = element_line(color = "lightgrey"),
          axis.ticks.y = element_blank(),
          plot.tag = element_blank(),      # an inset is never tagged
          plot.margin = margin(0, 0, 0, 0))

pB <- main + inset_element(inset, left = 0.48, bottom = 0.60,
                           right = 0.98, top = 0.98, align_to = "panel")

saveRDS(pB, file.path(OBJ_DIR, "B.rds"))
message("Panel B built.")
