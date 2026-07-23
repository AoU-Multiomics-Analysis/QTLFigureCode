# Panel F — median genes per individual by percent-change bin
source("fig2/00_setup.R")
median_pc <- readRDS(file.path(DATA_DIR, "F.rds"))

# Grey highlight boxes over the deletion-like / duplication-like bins.
DelBoxStart <- 0.5;  DelBoxEnd <- 2.45
DupBoxStart <- 7.5;  DupBoxEnd <- 13.45
DupBoxMid   <- DupBoxStart + (DupBoxEnd - DupBoxStart) / 2
max_median  <- max(median_pc$q75_genes)

pF <- median_pc %>%
  ggplot(aes(reorder(pct_bin, lower), median_genes)) +
    geom_rect(xmin = DelBoxStart, xmax = DelBoxEnd, ymin = -1, ymax = max_median + 100,
              inherit.aes = FALSE, fill = "lightgrey", alpha = 0.05) +
    geom_rect(xmin = DupBoxStart, xmax = DupBoxEnd, ymin = -1, ymax = max_median + 100,
              inherit.aes = FALSE, fill = "lightgrey", alpha = 0.05) +
    geom_col(width = 0.8, alpha = 0.6, linewidth = 0.4, color = "#462255FF", fill = "#462255FF") +
    geom_errorbar(aes(ymin = q25_genes, ymax = q75_genes), alpha = 0.7, width = 0.2) +
    facet_wrap(~gene_type, nrow = 2, strip.position = "right") +
    geom_text(x = 1.5,       y = max_median + 250, label = "Del-Like", size = ANN_MD,
              color = "grey35", fontface = "bold") +
    geom_text(x = DupBoxMid, y = max_median + 250, label = "Dup-Like", size = ANN_MD,
              color = "grey35", fontface = "bold") +
    geom_segment(aes(x = 0.6, xend = 2.25, y = 1000, yend = 1000),
                 inherit.aes = FALSE, color = "black", linewidth = 0.25) +
    geom_text(x = 1.425, y = 1200, label = "5 Genes", size = ANN_MD,
              color = "grey35", fontface = "bold") +
    geom_segment(aes(x = DupBoxStart + 0.5, xend = DupBoxEnd - 0.5, y = 1000, yend = 1000),
                 inherit.aes = FALSE, color = "black", linewidth = 0.25) +
    geom_text(x = DupBoxMid, y = 1200, label = "108 Genes", size = ANN_MD,
              color = "grey35", fontface = "bold") +
    labs(x = "Percent change bin", y = "Median genes per individual") +
    theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.x = element_text(angle = 90))

saveRDS(pF, file.path(OBJ_DIR, "F.rds"))
message("Panel F built.")
