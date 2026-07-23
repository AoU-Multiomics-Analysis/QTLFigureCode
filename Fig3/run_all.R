# =============================================================================
# Figure 3 — end-to-end runner. From the project root:  Rscript fig3/run_all.R
#   1. panels  : each stages the inputs it needs (fig3/inputs.tsv) and builds
#                the derived tables it needs (memoized in fig3/data/)
#                -> fig3/objects/<X>.rds
#   2. assemble: compose + save Figure3.pdf / Figure3.png
#
# To rebuild one panel, run just that panel script then fig3/assemble.R.
# The derived-data cache is keyed by NAME, not by content: clear fig3/data/
# after changing an input file.
#
# Supplementals are standalone:
#   Rscript fig3/supp_coloc_qc.R
#   Rscript fig3/supp_prop_gene_coloc.R
# =============================================================================
for (p in LETTERS[1:5]) source(sprintf("fig3/panels/panel%s.R", p))
source("fig3/assemble.R")
message("run_all complete.")
