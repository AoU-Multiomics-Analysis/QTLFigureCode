# Panel A — colocalizing-gene counts + credible-set coloc rate, by trait category
# Inputs: coloc_gene, coloc_rate
source("fig4/00_setup.R")

coloc_plot <- coloc_by_trait() %>%
  filter(trait_category != "neuropsychiatric", trait_category != "renal_urological") %>%
  filter(count > 5) %>%
  group_by(trait_category) %>%
  arrange(desc(count)) %>%
  filter(row_number() < 11) %>%
  build_plot() +
  theme(axis.text.x = element_text(size = 10))

saveRDS(coloc_plot, file.path(OBJ_DIR, "A.rds"))
message("Panel A built.")
