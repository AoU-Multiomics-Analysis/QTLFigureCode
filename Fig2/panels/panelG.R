# Panel G — validated dosage-like events (dumbbell)
source("fig2/00_setup.R")
outliers_df <- readRDS(file.path(DATA_DIR, "G.rds"))
outlier_pal2 <- c("Expression Underlier" = "#D44D5CFF", "Expression Overlier" = "#046E8FFF")

# Legend anchor, computed from the data. annotate() evaluates eagerly, so these
# are baked into the saved object -- assemble.R never needs them (see panelF.R
# for the lazy-aes() trap this avoids).
XMax     <- max(outliers_df$`Expression Overlier`)
LegX     <- XMax * 0.86
LegTextX <- XMax * 0.895

# Dumbbell drawn manually (segment + two end points) — no ggalt dependency.
pG <- ggplot(outliers_df, aes(y = label)) +
  geom_segment(aes(x = `Expression Underlier`, xend = `Expression Overlier`, yend = label),
               color = "grey65", linewidth = 2) +
  geom_point(aes(x = `Expression Underlier`), color = outlier_pal2["Expression Underlier"], size = 10) +
  geom_point(aes(x = `Expression Overlier`),  color = outlier_pal2["Expression Overlier"],  size = 10) +
  geom_text(aes(x = `Expression Underlier`, label = scales::comma(`Expression Underlier`)),
            hjust = 1.35, size = ANN_MD) +
  geom_text(aes(x = `Expression Overlier`, label = scales::comma(`Expression Overlier`)),
            hjust = -0.5, size = ANN_MD) +
  scale_x_continuous(labels = scales::comma, expand = expansion(mult = c(0.08, 0.34))) +
  # Legend sits in the bottom-right corner: y < 1 is below the lowest category,
  # so the band is empty across the full width. Anchored to the data range
  # rather than hard-coded x, so it stays in the corner if the counts change.
  annotate("point", x = LegX,      y = 0.75, color = "#D44D5CFF", size = 5) +
  annotate("text",  x = LegTextX,  y = 0.75, label = "Deletion-like",    hjust = 0, size = ANN_MD) +
  annotate("point", x = LegX,      y = 0.55, color = "#046E8FFF", size = 5) +
  annotate("text",  x = LegTextX,  y = 0.55, label = "Duplication-like", hjust = 0, size = ANN_MD) +
  labs(x = "Validated dosage-like events", y = NULL) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(color = "black", margin = margin(r = 25)),
        plot.margin = margin(t = 6, r = 10, b = 6, l = 4))

saveRDS(pG, file.path(OBJ_DIR, "G.rds"))
message("Panel G built.")
