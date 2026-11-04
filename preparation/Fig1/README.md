These scripts prepare the current Figure 1 and related supplement inputs. Run them from the repository root through the preparation runner. Supply the local input roots. All outputs go to `QTL_PREP_ROOT`.

| Step | Source and result |
| --- | --- |
| `fig1_gene_reference` | Read the GENCODE v48 GTF. Extract stable gene IDs, names, biotypes and coordinates. |
| `fig1_novelty` | Read raw COMB tensorQTL results and GTEx, INTERVAL, TOPMed and UKBB reference results. Rebuild the significant reference gene sets and current novelty calls. |
| `fig1_circos` | Prepare lead points, radial caps, selected labels, assay dots and axes from the current novel targets. |
| `fig1_categories` | Map NCG, ClinGen, gnomAD, TF, EpiFactors, RNA-processing, drug and GO annotations. Add HGNC histones and the membrane/signaling union. Count stable gene IDs once per assay and category. |
| `fig1_finemapping` | Aggregate raw variant rows into credible sets. Prepare set sizes, trait counts, lead PIPs, tied-lead counts, main D/E and S11 inputs. |
| `fig1_supplements` | Prepare S12 size summaries and gene-cluster standard errors. Prepare S13 bounded PIP densities and summary values. |
| `fig1_population_discovery` | Calculate the S4 BH FDR within each of the 21 permutation files. Count the common reference-cohort samples. |
| `fig1_cohort_counts` | Preserve supplied pre-QC aggregate counts. Count post-QC RNA samples and protein matrix/covariate header samples. |

Main Figure 1 novelty uses the source `qval < 0.05`. It excludes the paper's pseudogene, IG, TEC and TR biotypes, the artifact annotation, and pQTL lncRNAs. The reference comparison stays within the assay. The UKBB comparison uses protein targets from the ST9 cis-gene column. It expands combined targets and resolves identity with ST3 UniProt accessions, GENCODE and HGNC. It does not use the lead-variant gene annotation as the protein target.

S4 retains the original panel's separate rule: calculate BH FDR from `pval_beta` within each file. The common RNA cohort supplies its x-axis sample counts across assays. These are not assay-specific sample counts.

Set `QTL_FIG1_EAS_EQTL_FILE` to the current November 2 EAS eQTL download, or stage it at `Fig1Data/PopulationLevelFinemapping/eQTL/EAS_20261102.eQTL_SusieMerged.tsv` under the input root. The older `AoU_EAS_9k_eQTL_SusieMerged.tsv.gz` is not a fallback. The current download is gzip CSV despite its `.tsv` filename. The reader checks its compression signature.

The cohort step accepts `QTL_FIG1_PROTEIN_BED_FILE` and `QTL_FIG1_PROTEIN_COVARIATES_FILE`. Without these overrides, stage the two files at the `reference/` paths listed in the preparation manifest. The script writes aggregate counts only. The source participant availability export for the pre-QC table is unavailable; those supplied aggregate counts are retained.

The gene-reference and novelty steps can run without saved answer files. When these files are present, the scripts compare the prepared results with them. When they are absent, the novelty audit records `not_run` and null comparison values. Manual labels are inputs to the separate circos step.

The local checks matched 44,421 novelty calls, 13,741 circos points, all 21 category cells and their gene memberships, all 457,883 credible sets and lead PIPs, the S12/S13 summaries, and all 21 S4 summaries. The scripts do not rerun tensorQTL or SuSiE. Final SVG label spacing and icon artwork remain part of the archived figure composition.

The original script hashes are in the preparation manifest. No source tables, participant records or downloaded annotations belong in the repository. Run the synthetic fixture checks with `Rscript tests/test_preparation_fig1.R`.
