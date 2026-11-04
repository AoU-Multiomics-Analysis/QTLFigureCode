source("preparation/common.R")
source("preparation/Fig3/methods.R")
source("preparation/Fig3/symmetric_methods.R")
source("preparation/Fig3/current_methods.R")
source("preparation/Fig3/result_io.R")
prep_result_runtime("Fig3")

# Reconstruct both sides of each tested pair from source SuSiE tables. Coding
# calls require canonical CDS of the target gene, SpliceAI max <=0.2, and no
# cCRE overlap. Current panel E uses recoverable eQTL–pQTL colocalizing variants.
message("Prepare the symmetric credible-set background and annotation models")
enrichment <- coloc_maf_enrichment()
prep_write_tsv(as_tibble(enrichment), "source", "Fig3/data/current_panel_E.tsv")
message("Figure 3E data prepared; no sequence model or FastENLOC job was run")
