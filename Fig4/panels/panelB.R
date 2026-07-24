# Panel B — R^2 increase (%) contributed by rare variants
# Inputs: r2_twas
source("fig4/00_setup.R")

R2_combined <- read_any(input_path("r2_twas")) %>%
  mutate(H2_upper = (SE * 1.96) + H2, adj_R2 = r.squared / H2) %>%
  filter(Pval < 1e-5) %>%
  select(phenotype_id, AF_threshold, adj_R2, H2, Pval) %>%
  pivot_wider(names_from = AF_threshold, values_from = adj_R2) %>%
  mutate(difference = (`0` - `0.01`))

R2_dist <- R2_combined %>%
  ggplot(aes(x = 100 * difference)) +
  geom_histogram(color = "white") +
  ylab("Number of genes") +
  xlim(-20, 40) +
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14)) +
  labs(x = expression(R^2 ~ increase ~ "(%)"))

saveRDS(R2_dist, file.path(OBJ_DIR, "B.rds"))
message("Panel B built.")
