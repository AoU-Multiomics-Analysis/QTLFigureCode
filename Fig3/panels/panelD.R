# Panel D — eQTL vs pQTL posterior effect size, faceted by MAF bin.
# Inputs: enloc_variant, annotations, susie_eqtl, susie_pqtl
source("fig3/00_setup.R")

binned <- coloc_model_input()[["eQTL-pQTL"]] %>%
  filter(!is.na(max_af)) %>%
  mutate(maf_bin = case_when(max_af >= .05 ~ "Common",
                             max_af >= .01 ~ "Low-Freq",
                             TRUE          ~ "Rare"),
         posterior_mean.x = log1p(posterior_mean.x),
         posterior_mean.y = log1p(posterior_mean.y))

lims <- quantile(abs(c(binned$posterior_mean.x, binned$posterior_mean.y)),
                 0.97, na.rm = TRUE)

slopes <- binned %>%
  group_by(maf_bin) %>%
  summarise(n     = dplyr::n(),
            slope = coef(lm(posterior_mean.y ~ posterior_mean.x))[["posterior_mean.x"]],
            .groups = "drop") %>%
  mutate(label = paste0("β = ", round(slope, 2), "\nN = ", scales::comma(n)))

# Every row here is already Coloc == TRUE (coloc_model_input filters), so the
# original's redundant filter(Coloc == TRUE) inside geom_point is gone.
pD <- ggplot(binned, aes(x = posterior_mean.x, y = posterior_mean.y)) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "grey85") +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "grey85") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              linewidth = 0.35, color = "grey65") +
  geom_point(alpha = 0.6, size = 4, color = "#6EC5C3") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.7,
              color = "red", fill = "red", alpha = 0.12) +
  geom_label(data = slopes, aes(x = -Inf, y = Inf, label = label),
             inherit.aes = FALSE, hjust = -0.05, vjust = 1.15,
             size = SZ_ANNOT, label.size = 0, fill = "white", alpha = 0.85) +
  facet_wrap(~maf_bin, nrow = 1) +
  coord_cartesian(xlim = c(-lims, lims), ylim = c(-lims, lims)) +
  labs(x = "Expression Posterior", y = "Protein Posterior") +
  theme(panel.spacing.x = unit(0.8, "lines"),
        plot.margin = margin(5.5, 8, 5.5, 5.5))

saveRDS(pD, file.path(OBJ_DIR, "D.rds"))
message("Panel D built.")
