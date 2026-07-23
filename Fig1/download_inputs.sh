#!/usr/bin/env bash
# =============================================================================
# Provenance record for Figure 1 inputs — NOT executed automatically.
# -----------------------------------------------------------------------------
# These are the Google Cloud Storage (and local Terra) origins of every raw
# input consumed by fig1/01_prep_data.R, captured from the original notebook
# (QTLAnalysisUtils/notebooks/QTL_paper_Fig1.ipynb).
#
# The data is All of Us CONTROLLED-ACCESS; downloading requires the appropriate
# Researcher Workbench credentials and `gsutil`. To stage inputs, uncomment the
# relevant lines (or copy the URIs) and run from the project root. The commands
# below place files into the input_data/ layout expected by 01_prep_data.R.
# =============================================================================
set -euo pipefail
BUCKET="gs://fc-secure-b8771cfd-5455-4292-a720-8533eb501a93"

# mkdir -p input_data/{reference,external,annotated}
# mkdir -p input_data/permutations/{eqtl,sqtl,pqtl}
# mkdir -p input_data/finemapping/{eqtl,sqtl,pqtl}

# --- reference/ --------------------------------------------------------------
# gsutil cp "$BUCKET/eQTL_9k_set/sample_lists/COMB_eQTL_9k_sample_list.tsv" input_data/reference/
# Local (Terra) inputs — not in the bucket above; stage from their source:
#   input_data/reference/gencode.v48.annotation.gtf.gz   (/home/jupyter/multiomics_qtls/edit/)
#   input_data/reference/echo_v4_r2.ancestry_preds.tsv   (/home/jupyter/multiomics_qtls/edit/)

# --- permutations/eqtl/ ({ANCESTRY}.cis_qtl.txt.gz) --------------------------
# submissions/1efb1546-1a5a-4b3d-a4a0-9f36a111d8b2/tensorqtl_cis_permutations_workflow/
#   AFR: .../e15f1806-74dd-4424-96d8-84dfb4014b8b/call-tensorqtl_cis_permutations/AFR.cis_qtl.txt.gz
#   AMR: .../3c72146a-a8eb-4086-8343-261def644141/call-tensorqtl_cis_permutations/AMR.cis_qtl.txt.gz
#   COMB:.../5d07a1a2-f3ea-4781-8fdf-f15cf5102ba6/call-tensorqtl_cis_permutations/COMB.cis_qtl.txt.gz
#   EAS: .../4bbcdb39-8c87-424d-bc3e-23fedea7413b/call-tensorqtl_cis_permutations/EAS.cis_qtl.txt.gz
#   EUR: .../f87a235b-5a3a-43b6-8fb4-868beb2a48ba/call-tensorqtl_cis_permutations/EUR.cis_qtl.txt.gz
#   MID: .../31ffd78f-f879-42a8-9dd8-7ca3f5914585/call-tensorqtl_cis_permutations/MID.cis_qtl.txt.gz
#   SAS: .../624546ec-8dac-453d-9de9-2b531705e7b2/call-tensorqtl_cis_permutations/SAS.cis_qtl.txt.gz

# --- permutations/pqtl/ (AoU_10750_{ANCESTRY}_pQTL.cis_qtl.txt.gz) -----------
# submissions/ff89514b-ca9f-4306-8f0b-2b3fe2784c6b/... for AFR/AMR/EAS/EUR/MID/SAS
#   COMB is under submissions/7b9df902-758b-474a-8a62-1354eacbed9d/...
#   (call-tensorqtl_cis_permutations/AoU_10750_{ANCESTRY}_pQTL.cis_qtl.txt.gz)

# --- permutations/sqtl/ ({ANCESTRY}.cis_qtl.txt.gz) --------------------------
# submissions/3706a9f9-74c5-4662-aae2-b71875852223/tensorqtl_cis_permutations_workflow/
#   AFR/AMR/COMB/EAS/EUR/MID/SAS under the per-call subfolders (see notebook).

# --- finemapping/eqtl/ (AoU_{ANCESTRY}_9k_eQTL_SusieMerged.tsv.gz) -----------
# AggregateSusieWorkflow/call-AggregateSusie/ for EAS/AMR/AFR/EUR/SAS/COMB/MID
#   (submission IDs per ancestry — see notebook cell "Load fine-mapping data").

# --- finemapping/sqtl/ (AoU_{ANCESTRY}_9k_sQTL_SusieMerged.annotated.tsv) ----
# AggregateSusieWorkflow/call-AnnotateSusie/ for COMB/SAS/EUR/EAS/AMR/AFR/MID.

# --- finemapping/pqtl/ -------------------------------------------------------
#   input_data/finemapping/pqtl/pQTL_fm_data_10k.tsv   (local Terra working file)

# --- annotated/ (COMB tables for functional-enrichment, Panel H) -------------
# gsutil cp "$BUCKET/submissions/8e000b9e-17db-44f1-b7eb-4a14db8bd384/AnnotateSusie/76fc2dc3-f7b8-4c66-b0c4-645eac86caeb/call-AnnotateSusieTask/AoU_COMB_9k_eQTL_SusieMerged.annotated.tsv" input_data/annotated/
# gsutil cp "$BUCKET/submissions/77cfa98c-8dbf-4d38-9e4b-9e9b23d2308f/AggregateSusieWorkflow/93c8280a-f41d-4584-b25c-0ad9ca3910a9/call-AnnotateSusie/AoU_COMB_9k_sQTL_SusieMerged.annotated.tsv" input_data/annotated/
#   input_data/annotated/pQTL_susie_annotated.tsv       (local Terra working file)

# --- external/ (reference consortia) ----------------------------------------
# gsutil cp "$BUCKET/INTERVAL_eQTL_summary.xlsx"                 input_data/external/
# gsutil cp "$BUCKET/INTERVAL_cis_sQTL_summary.xlsx"            input_data/external/
# gsutil cp "$BUCKET/Whole_blood.maf0001.100.permutations.txt"  input_data/external/   # TOPMed eQTL
# Local (Terra) GTEx v10 Whole Blood:
#   input_data/external/Whole_Blood.v10.eGenes.txt.gz  (/home/jupyter/GTEx_Analysis_v10_eQTL_updated/)
#   input_data/external/Whole_Blood.v10.sGenes.txt.gz  (/home/jupyter/GTEx_Analysis_v10_sQTL_updated/)

# --- Deferred (see NOTES.md): UKBB pQTL reference for the novel-pGene comparison
#   input_data/external/<UKBB pQTL supplementary table>.xlsx    (path TBD)
