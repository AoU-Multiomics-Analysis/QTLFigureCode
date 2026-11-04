# %% Source-result readers and transformations
#| Start with GWAS colocalization, TWAS, ClinGen and DGIdb result tables. These calculations run in this R session.
# @include preparation/Fig4/methods.R
# @include preparation/Fig4/current_methods.R

# @include notebooks/current/plots/Fig4.R

# %% Locate the Figure 4 source results
#| Resolve the original result-table manifest. Session memory is the only cache.
prep_result_runtime("Fig4")
message("Figure 4 source-result readers are ready.")

# %% Count colocalizing genes and join union rates
#| Count distinct gene–trait pairs with any_coloc. Join the lenient union rate, remove the two excluded categories, require more than five genes, and keep at most ten traits per category. The full harmonized gene table is required.
fig4_gene_relative <- "Fig4Data/high_level_coloc_outputs/harmonized/harmonized_coloc.gene.tsv.gz"
fig4_gene_available <- file.exists(file.path(prep_root("input"), fig4_gene_relative))
a <- NULL
if (fig4_gene_available) {
  coloc_gene <- read_any(input_path("coloc_gene"))
  coloc_rates <- read_any(input_path("coloc_rate"))
  coloc_traits <- build_plot_data(coloc_gene, coloc_rates, min_coloc_genes = 1, rate_type = "lenient")
  coloc_display <- fig4_current_trait_data(coloc_traits)
  a <- fig4_plot_traits(coloc_display)
  print(a)
} else {
  message("Figure 4A skipped: missing full harmonized gene result: ", fig4_gene_relative,
          ". The source code is shown above. Panels B–D will be calculated.")
}

# %% Calculate rare-variant prediction contrasts
#| Retain genes with heritability P < 1e-5. Divide prediction R-squared by heritability, then subtract the AF >= 0.01 prediction from the AF >= 0 prediction.
r2_source <- read_any(input_path("r2_twas"))
r2_combined <- fig4_r2_data(r2_source)
b <- fig4_plot_r2(r2_combined)
print(b)

# %% Select conditional TWAS gene–trait pairs
#| Remove failed results and duplicate molecular-trait rows. Exclude rows with abs(z_rare_cond)/max(abs(rare_zscore), 0.5) > 3 and abs(Vcor) > 0.3. For each gene–trait pair, select the greatest GLCP, then the smallest conditional rare-variant P value. Require finite matching covariates.
conditional_twas <- read_any(input_path("conditional_twas"))
qq_match_input <- prepare_qq_match_source(conditional_twas)
message("Prepared conditional TWAS matching covariates.")

# %% Match controls and calculate QQ coordinates
#| Use nearest-neighbor propensity matching with seed 20260713, exact trait matching, ratio 1, no replacement, and a 0.2 standardized caliper. Use log TPM, log shet, log rare-variant count and abs(Vcor). Rank P values within each test and group. Cap display values at 40 for common variants and 12.5 for rare variants.
qq_match_result <- match_colocalized_pairs(qq_match_input, seed = 20260713, ratio = 1,
                                           caliper = 0.2, replace = FALSE)
qq_matched <- qq_match_result$data
qq_blood <- calculate_qq_data(qq_matched, traits = BLOOD_TRAITS)
qq_display <- qq_blood %>% mutate(group = case_when(str_detect(group, "Matched") ~ "Null", TRUE ~ group))
c <- fig4_plot_qq(qq_display)
print(c)

# %% Define ClinGen and drug-target gene sets
#| Keep definitive or strong ClinGen evidence. Apply the visible immune/hematological disease-name rules. Use upper-case DGIdb gene symbols for druggable and approved target sets.
clingen_source <- arrow::read_parquet(input_path("clingen")) %>%
  select(targetFromSourceId, allelicRequirements, confidence, score, studyId, targetId, diseaseFromSource)
drug_source <- read_any(input_path("drug_interactions"))
gs <- derive_gene_sets(clingen_source, drug_source)
message("Prepared ClinGen and DGIdb gene sets.")

# %% Adjust TWAS P values within each source path
#| Retain colocalizing results after duplicate and artifact removal. Apply Benjamini–Hochberg adjustment separately within each TWAS_path. A rare result requires conditional FDR < 0.10 and unadjusted rare P < 0.05.
coloc_adjusted <- conditional_twas %>%
  distinct(gene, coloc_trait, molecular_trait_id, .keep_all = TRUE) %>%
  mutate(inflation = abs(z_rare_cond) / pmax(abs(rare_zscore), 0.5),
         artifact_primary = coalesce(inflation > 3 & abs(Vcor) > 0.3, FALSE)) %>%
  filter(!artifact_primary, any_coloc == TRUE) %>%
  group_by(TWAS_path) %>%
  mutate(fdr_rare_cond = p.adjust(p_rare_cond, method = "BH"),
         fdr_common = p.adjust(common_pvalue, method = "BH"),
         rare_zscore_fdr = p.adjust(rare_pvalue, method = "BH"),
         fdr_full = p.adjust(full_pvalue, method = "BH")) %>%
  ungroup() %>%
  mutate(genuine_rare = p_rare_cond < 0.05 & rare_pvalue < 0.05,
         t_fdrcond = fdr_rare_cond < 0.10 & rare_pvalue < 0.05,
         t_fdrboth = fdr_rare_cond < 0.10 & rare_zscore_fdr < 0.10,
         gene = factor(gene), coloc_trait = factor(coloc_trait),
         druggable = toupper(as.character(Description)) %in% gs$druggable,
         approved = toupper(as.character(Description)) %in% gs$approved)
stopifnot(nrow(filter(coloc_adjusted, abs(Vcor) > 0.3, inflation > 3)) == 0L)
message("Applied TWAS-path FDR adjustment after artifact removal.")

# %% Summarize gene-level TWAS evidence
#| Use the maximum absolute rare and common Z scores per gene, mean TPM, the first shet value, and whether the gene has any qualifying rare or common result. Add ClinGen membership and retain nonmissing log10(shet).
gtd <- coloc_adjusted %>%
  group_by(gene, Description) %>%
  summarise(max_rare = max(abs(z_rare_cond), na.rm = TRUE),
            max_common = max(abs(common_zscore), na.rm = TRUE),
            TPM = mean(mean_TPM, na.rm = TRUE), druggable = first(druggable),
            approved = first(approved), ever_common = any(fdr_common < 0.1, na.rm = TRUE),
            ever_rare = any(t_fdrcond), shet = first(shet), .groups = "drop") %>%
  mutate(log_shet = log10(shet), log_tpm = log10(TPM + 1),
         MendelianGene = as.character(gene) %in% gs$clingen_cleaned$targetId,
         HemeImmuneGene = as.character(gene) %in% gs$heme_immune) %>%
  filter(!is.na(log_shet))
message("Prepared gene-level TWAS evidence for enrichment.")

# %% Fit gene-set enrichment models
#| Fit binomial models for each gene set with continuous and binary TWAS predictors, log10(shet), and log2(TPM). Exponentiate coefficients and profile confidence limits. Plot only the continuous rare and common Z-score terms.
enrichment <- fig4_enrichment_data(gtd)
enrichment_results <- enrichment$all_fits
enrichment_display <- enrichment$plot_data
d <- fig4_plot_enrichment(enrichment_display)
print(d)

# %% Render Figure 4
#| Assemble the calculated panels at 184 x 190 mm. If the full harmonized gene input is absent, panel A remains empty.
fig <- fig4_assemble(a, b, c, d)
save_plot(fig, "Fig4/Figure4_rebuilt", 184, 190)
print(fig)

# %% Reference figure
#| Complete figure layout.
show_reference(4)
