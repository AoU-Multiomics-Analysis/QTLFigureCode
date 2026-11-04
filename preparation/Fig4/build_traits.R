source("preparation/common.R")
source("preparation/Fig4/methods.R")
source("preparation/Fig4/current_methods.R")
source("preparation/Fig3/result_io.R")
prep_result_runtime("Fig4")

# Requires the harmonized gene result with trait_category and any_coloc. A
# final membership table without these fields cannot replace that source.
message("Prepare trait-level gene counts and lenient union rates")
traits <- coloc_by_trait()
prep_write_tsv(fig4_current_trait_data(traits), "source", "Fig4/data/current_panel_A.tsv")
