# Panel A (labelled variant) — same data and layout as panelA.R, plus a short
# label on every trait, right-aligned into the dead space to the right of the
# bars. Smaller than build_plot()'s default label size: 70 rows share the panel
# height, so the type has to fit a ~9pt row pitch. Feeds Figure4_labeled.pdf;
# Figure4.pdf still uses the plain A.rds.
# Inputs: coloc_gene, coloc_rate
source("fig4/00_setup.R")

coloc_plot <- coloc_by_trait() %>%
  filter(trait_category != "neuropsychiatric", trait_category != "renal_urological") %>%
  filter(count > 5) %>%
  group_by(trait_category) %>%
  arrange(desc(count)) %>%
  filter(row_number() < 11) %>%
  build_plot(label_top = Inf, label_size = 2.2) +
  theme(axis.text.x = element_text(size = 10))

saveRDS(coloc_plot, file.path(OBJ_DIR, "A_labeled.rds"))
message("Panel A (labelled) built.")
