# %% Read and label the original FastENLOC and SuSiE results
#| Read the original gene, credible-set, and variant result tables. The functions below apply Bayesian FDR 5% separately to each assay pair. SuSiE tables supply the shared-gene set, allele frequencies, and posterior effects. These steps process completed results; they do not rerun FastENLOC or fine-mapping.
# @include preparation/Fig3/methods.R
prep_result_runtime("Fig3")
message("Figure 3 source readers and FDR rules are ready.")

# %% Define the symmetric annotation background
#| FastENLOC lists one side of a credible-set pairing. Recover the other side from the corresponding SuSiE table. Coding calls require canonical CDS for the target gene, SpliceAI maximum at most 0.2, and no cCRE overlap. The GTF reader in the next cell replaces the platform-specific reader defined here.
# @include preparation/Fig3/symmetric_methods.R
message("Symmetric credible-set and annotation methods are ready.")

# %% Define the current summaries and model formulas
#| Fold each ancestry frequency to MAF before taking the maximum across ancestries. Keep missing frequencies. The effect models retain the current signed eQTL-pQTL formula and the absolute-effect formulas for the other pairs. Standardize negative log10 MAF for the interaction models.
# @include preparation/Fig3/current_methods.R
message("Current Figure 3 summaries and model formulas are ready.")

# %% Define plots from the calculated frames
#| Construct all five plots from the new frames. Plot builders do not read saved data or plot objects.
# @include notebooks/current/plots/Fig3.R
message("Figure 3 plot builders are ready.")

# %% Gene-level sharing
#| Use genes fine-mapped in all three assays as the denominator. For each gene, count the assay-pair FDR calls, then calculate the fraction of genes in each count. The current inputs contain 755 shared genes.
fm <- susie_finemapping()
gene <- coloc_gene_level()
layer_counts <- fig3_gene_counts(gene, fm$genes_all_three_assays)
stopifnot(sum(layer_counts$n) == length(fm$genes_all_three_assays),
          abs(sum(layer_counts$prop) - 1) < 1e-12)
panels <- list(A = fig3_panel_a(layer_counts))
save_plot(panels$A, "Fig3/panel_A", 130, 90)
print(panels$A)

# %% Colocalized signal components
#| Build an undirected graph from FDR-positive credible-set pairs. Count its connected components per gene. The main plot uses all genes with a component. The inset restricts the single-assay credible-set denominator to the shared-gene set. These are different populations.
cs <- coloc_cs_level()
components <- coloc_cs_components()
signal_summaries <- fig3_signal_counts(components, fm$genes_all_three_assays)
panels$B <- fig3_panel_b(signal_summaries)
save_plot(panels$B, "Fig3/panel_B", 150, 100)
print(panels$B)

# %% Join posterior effects and refit interaction models
#| Join variant annotations by variant and gene-relative annotations by variant plus gene. For each variant-gene key, use the largest absolute posterior effect within each assay. Retain FDR-positive variant calls, fit the three current linear models, and plot their interaction coefficients with 95% confidence limits.
variant <- coloc_variant_level()
models <- coloc_model_input()
interaction_results <- fig3_interaction_data(models)
panels$C <- fig3_panel_c(interaction_results)
stopifnot(nrow(panels$C$data) == 3L)
save_plot(panels$C, "Fig3/panel_C", 140, 100)
print(panels$C)

# %% Transform expression and protein effects by MAF
#| Apply log1p to the signed eQTL-pQTL effects. Use Rare below 0.01, Low-Freq from 0.01 to below 0.05, and Common from 0.05. Fit a linear slope within each group. Limit the display to the 97th percentile of absolute transformed effects; fit the models before this display limit. The current inputs contain 1,722 effect rows.
effects <- fig3_effect_data(models)
panels$D <- fig3_panel_d(effects)
save_plot(panels$D, "Fig3/panel_D", 180, 100)
print(panels$D)

# %% Refit variant annotation models
#| Use recoverable colocalizing eQTL-pQTL variants from both assay sides. Collapse repeated variants with geometric-mean PIP weights and combined assay membership. For each annotation, fit log10 MAF with PIP weights and adjust for assay membership and PIP. Negate the annotation coefficient and use base 10 for the MAF fold change. Apply BH adjustment across fitted annotations; show adjusted P below 0.05 and at least five annotated variants.
enrichment <- coloc_maf_enrichment()
panels$E <- fig3_panel_e(enrichment)
save_plot(panels$E, "Fig3/panel_E", 150, 100)
print(panels$E)

# %% Assemble the five calculated panels
#| Keep one tag for panel B and its inset. Use the current panel order and layout.
fig3 <- fig3_assemble(panels)
save_plot(fig3, "Fig3/Figure3_rebuilt", 184, 130)
print(fig3)

# %% Reference figure
#| Complete figure layout.
show_reference(3)
