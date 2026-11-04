# Pure fixtures for the current Figure 1 preparation rules. No source records.
suppressPackageStartupMessages(library(tidyverse))
source("tests/helpers.R")
source("preparation/Fig1/functions.R")

genes <- tibble(gene_id = paste0("ENSG", 1:6), gene_name = LETTERS[1:6],
                gene_type = c("protein_coding", "lncRNA", "artifact", "processed_pseudogene", "snRNA", "protein_coding"))
raw <- tibble(assay = c("eQTL", "pQTL", "eQTL", "eQTL", "sQTL", "eQTL"),
              phenotype_id = paste0(genes$gene_id, ".1"), variant_id = "chr1:10_A_G",
              source_qval = c("1e-310", "0.001", "0.001", "0.001", "0.049", "0.05"),
              source_pval_nominal = "0", source_pval_beta = "0.001")
known <- tibble(assay = "eQTL", gene_id = "ENSG5", study = "GTEx")
novel <- fig1_classify_novelty(raw, genes, known, measured_ids = character())
ok("use source q values; exclude artifact, pseudogene, protein lncRNA and q=0.05",
   identical(novel$novel_for_figure, c(TRUE, FALSE, FALSE, FALSE, TRUE, FALSE)))
ok("reference comparison stays within the assay", novel$novel_for_figure[[5]])
ok("keep source zero and subnormal P-value text", novel$source_qval[[1]] == "1e-310" &&
     novel$source_pval_nominal[[1]] == "0")
expect_error("reject duplicate source trait keys", fig1_classify_novelty(bind_rows(raw, raw[1, ]), genes, known))

hgnc <- tibble(status = "Approved", symbol = c("NEW", "OTHER", "ALIAS1", "ALIAS2"),
               ensembl_gene_id = c("ENSG1", "ENSG2", "ENSG3", "ENSG4"),
               prev_symbol = c("OLD", NA, NA, NA), alias_symbol = c(NA, NA, "AMB", "AMB"),
               uniprot_ids = c("P1", "P2", NA, NA))
protein_genes <- tibble(gene_id = c("ENSG1", "ENSG2", "ENSG7"),
                        gene_name = c("NEW", "OTHER", "NEW"), gene_type = "protein_coding")
panel <- tibble(assay_target = c("old assay", "shared symbol"), gene_name = c("OLD", "NEW"),
                uniprot_id = c("P2", NA))
mapping <- fig1_ukbb_target_map(protein_genes, hgnc, panel)
ok("UniProt resolves a conflicting historical protein identity",
   identical(mapping$map$gene_id[mapping$map$gene_name == "OLD"], "ENSG2"))
ok("preserve verified overlapping GENCODE protein IDs",
   setequal(mapping$map$gene_id[mapping$map$gene_name == "NEW"], c("ENSG1", "ENSG7")))
ok("exclude an ambiguous HGNC alias", !"AMB" %in% mapping$map$gene_name)
st9 <- tibble(`cis/trans` = c("cis", "trans"), `Assay Target` = "assay",
               `cis gene` = c("OLD, NEW", "OTHER"), `Ensembl gene ID` = "ENSG999")
targets <- fig1_ukbb_cis_targets(st9, mapping$map)
ok("expand ST9 targets and never use the variant gene column",
   setequal(targets$mapped$gene_id, c("ENSG1", "ENSG2", "ENSG7")) &&
     !"ENSG999" %in% targets$mapped$gene_id)

variants <- tibble(molecular_trait_id = c("ENSG1.1", "ENSG1.1", "ENSG1.1", "ENSG1.1"),
                   cs_id = c("ENSG1.1_L1", "ENSG1.1_L1", "ENSG1.1_L2", "ENSG1.1_L2"),
                   cs_size = c(2L, 2L, 2L, 2L), variant = c("a", "b", "c", "d"),
                   pip = c(.9, .9, .95, .1))
sets <- fig1_credible_sets(variants, "eQTL", "COMB")
ok("one lead per set; report tied lead variants", nrow(sets) == 2 && sets$n_tied_lead_variants[[1]] == 2)
summary <- fig1_finemapping_summaries(sets)
ok("strict high-PIP cutoff and per-trait credible-set counts",
   summary$lead_summary$n_pip_above_09 == 1L && summary$traits$number_credible_sets == 2L)
expect_error("reject inconsistent credible-set size", fig1_credible_sets(
  variants %>% mutate(cs_size = c(2L, 3L, 2L, 2L)), "eQTL", "COMB"))
expect_error("reject an unresolved cloud path", fig1_assert_local_input("gs://bucket/file.tsv"))

membership <- tibble(assay = c("eQTL", "eQTL", "sQTL"), gene_id = c("ENSG1", "ENSG1", "ENSG1"),
                     gene_name = "A", cancer_driver = TRUE, lof_constrained = NA)
category <- fig1_category_counts(membership, c("cancer_driver", "lof_constrained"))
ok("category counts use distinct genes and preserve unknown membership",
   all(category$counts$n_genes[category$counts$category == "cancer_driver"] == 1L) &&
     all(category$counts$n_genes[category$counts$category == "lof_constrained"] == 0L))
nodes <- tibble(go_id = c("GO:1", "GO:2", "GO:3"))
edges <- tibble(child = c("GO:2", "GO:3"), parent = "GO:1", predicate = c("is_a", "regulates"))
ok("GO descendants use only is_a and part_of edges",
   setequal(fig1_go_descendants(nodes, edges, "GO:1"), c("GO:1", "GO:2")))
message("Current Figure 1 preparation fixtures passed.")
