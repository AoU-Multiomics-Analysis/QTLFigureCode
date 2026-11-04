source("preparation/common.R")
source("preparation/Fig4/methods.R")
source("preparation/Fig4/current_methods.R")
source("preparation/Fig3/result_io.R")
prep_result_runtime("Fig4")

message("Prepare rare-variant prediction contrasts")
prep_save_rds(fig4_r2_data(read_any(input_path("r2_twas"))), "source", "Fig4/data/current_panel_B.rds")

# The source result-table boundary already contains conditional-TWAS tests.
# Artifact removal and propensity matching are repeated with the original
# fixed seed, exact trait matching, covariates, ratio, and caliper.
message("Remove conditional-TWAS artifacts and match control gene–trait pairs")
matched <- qq_match_source()
qq <- calculate_qq_data(matched, traits = BLOOD_TRAITS)
prep_save_rds(qq, "source", "Fig4/data/current_panel_C.rds")

# Refit downstream gene-set enrichment after the original per-TWAS-path BH
# adjustments, gene-level collapse, and ClinGen/DGIdb joins.
message("Prepare gene sets and refit downstream enrichment models")
enrichment <- fig4_enrichment_data(gene_twas_data())
prep_save_rds(enrichment, "source", "Fig4/data/current_panel_D.rds")
prep_write_tsv(enrichment$plot_data, "source", "Fig4/data/current_panel_D.tsv")
message("Figure 4 B–D plotting data prepared; upstream TWAS was not rerun")
