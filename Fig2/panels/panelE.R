# Panel E — QTL-burden outlier enrichment
source("fig2/00_setup.R")
OutlierEnrichments <- readRDS(file.path(DATA_DIR, "E.rds"))
outlier_pal <- c("#D44D5CFF", "#046E8FFF")

pE <- OutlierEnrichments %>%
  filter(focal_lower_bound != -10) %>%
  ggplot(aes(focal_lower_bound, odds_ratio, color = type, fill = type, group = type)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_line(show.legend = FALSE) +
    geom_linerange(aes(ymin = ci_low, ymax = ci_high), show.legend = FALSE) +
    geom_point(size = 6, shape = 21, alpha = 0.7, show.legend = FALSE) +
    scale_fill_manual(values = outlier_pal) +
    scale_color_manual(values = outlier_pal) +
    scale_y_continuous(trans = "log2", breaks = 2^c(-4, 0, 4, 8),
                       labels = function(x) parse(text = paste0("2^", log2(x)))) +
    annotate("text", x = -100, y = 7, label = "Underlier", color = "#D44D5CFF", size = ANN_SM) +
    annotate("text", x =  175, y = 4, label = "Overlier",  color = "#046E8FFF", size = ANN_SM) +
    labs(y = "Outlier Enrichment") +
    theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.x = element_blank(), axis.title.x = element_blank())

saveRDS(pE, file.path(OBJ_DIR, "E.rds"))
message("Panel E built.")
