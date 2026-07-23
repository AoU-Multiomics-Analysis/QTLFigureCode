# Panel E — CDF of credible sets per eGene (COMB)
# Inputs: finemap_eqtl (all ancestries)
source("fig1/00_setup.R")

eQTL_finemapping <- load_input("finemap_eqtl")
eQTL_cs_number <- eQTL_finemapping %>%
  distinct(molecular_trait_id, cs_id, group, Assay)
ccdf <- ComputeCredibleSetCDF(eQTL_cs_number)

pE <- ccdf %>%
  filter(number_cs > 0, group == "COMB", number_cs <= 5) %>%
  ggplot(aes(y = prop_ge, x = number_cs)) +
    geom_step(color = "#CF597E", lwd = 3.5) +
    geom_point(shape = 21, size = 10, fill = "#CF597E") +
    scale_x_continuous(breaks = 1:5, labels = c("1", "2", "3", "4", "5+")) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    theme(
      plot.margin = margin(b = 20, unit = "pt"),
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    labs(x = "Number of credible sets",
         y = "eGenes with ≥ k \ncredible sets (%)")

saveRDS(pE, file.path(OBJ_DIR, "E.rds"))
message("Panel E built.")
