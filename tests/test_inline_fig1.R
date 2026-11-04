# Removing the visible builders must fail this test before any source data is read.
source("tests/helpers.R")
suppressPackageStartupMessages(library(tidyverse))
helper <- "notebooks/current/analysis/Fig1_helpers.R"
if (file.exists(helper)) source(helper)
ok("the visible Figure 1 builders exist", exists("fig1_size_contrasts", mode = "function"))
# AFR and COMB share genes. Their covariance must change the contrast SE.
sets <- tribble(~assay, ~population, ~cs_id, ~cs_size,
  "eQTL", "AFR", "ENSG1_L1", 1,
  "eQTL", "AFR", "ENSG2_L1", 3,
  "eQTL", "COMB", "ENSG1_L1", 2,
  "eQTL", "COMB", "ENSG2_L1", 6)
contrast <- fig1_size_contrasts(sets, "eQTL", c("COMB", "AFR"))
ok("shared genes contribute to the contrast uncertainty",
  isTRUE(all.equal(unname(contrast$standard_error[contrast$population == "COMB"]), 1)))
ok("the reference contrast has zero uncertainty",
  contrast$difference[contrast$population == "AFR"] == 0 &&
    contrast$standard_error[contrast$population == "AFR"] == 0)
# A row per variant must reduce to one record per set, with tied leads retained.
variants <- tribble(~molecular_trait_id, ~cs_id, ~cs_size, ~variant, ~pip,
 "ENSG1", "ENSG1_L1", 3L, "v1", .9,
 "ENSG1", "ENSG1_L1", 3L, "v2", .9,
 "ENSG1", "ENSG1_L1", 3L, "v3", .1)
cs <- fig1_credible_sets(variants, "eQTL", "COMB")
ok("one credible set retains its maximum PIP and tied leads",
  nrow(cs) == 1L && cs$lead_pip == .9 && cs$n_tied_lead_variants == 2L)
ok("the 0.9 threshold is strict", fig1_finemapping_summaries(cs)$lead_summary$n_pip_above_09 == 0L)
nodes <- tibble(go_id = c("GO:1", "GO:2", "GO:3", "GO:4"))
edges <- tribble(~child, ~parent, ~predicate,
 "GO:2", "GO:1", "is_a", "GO:3", "GO:2", "http://purl.obolibrary.org/obo/BFO_0000050",
 "GO:4", "GO:1", "regulates", "GO:1", "GO:3", "is_a")
ok("GO unions include is_a and part_of descendants and terminate at cycles",
 setequal(fig1_go_descendants(nodes, edges, "GO:1"), c("GO:1", "GO:2", "GO:3")))
