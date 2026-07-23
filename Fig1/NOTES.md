# Figure 1 — deferred analyses & notes

## Deferred: novel pGene comparison (Panel D)

Panel D currently shows **novel eGene + sGene** counts relative to external
consortia. The third bar — **novel pGenes** — is deferred because it depends on a
UKBB pQTL reference table that was never wired up in the original notebook
(`UKBB_pqtl_path` was referenced but never assigned; it is a latent bug rather
than a deliberate omission).

**To re-enable it** once the reference is available:

1. Stage the UKBB pQTL supplementary table into `input_data/external/`
   (e.g. `input_data/external/UKBB_pQTL_supplementary.xlsx`).
2. In `fig1/01_prep_data.R`, load it and compute `novel_pQTLs`:
   ```r
   UKBB_pQTLs <- LoadUKBBpQTL(in_path("external", "UKBB_pQTL_supplementary.xlsx"))
   novel_pQTLs <- QTL_data %>%
     filter(group == "COMB" & Assay == "Protein") %>%
     filter(qval < .05) %>%
     mutate(phenotype_id = str_remove_all(phenotype_id, ".*_|\\..*")) %>%
     left_join(gene_types, by = c("phenotype_id" = "gene_id")) %>%
     ExtractNovelQTLs(UKBB_pQTLs$gene_id) %>%
     filter(group != "lncRNA") %>%
     mutate(type = "pGene")
   ```
3. `bind_rows(novel_eQTLs, novel_sQTLs, novel_pQTLs)` and add `"pGene"` to the
   `type` factor levels in `novel_qtl_count`.

`LoadUKBBpQTL()` is already folded into `fig1/00_setup.R` (bug-fixed), so only the
input file and the prep block above are needed.

## MetaSoft supplemental

`fig1/supp_metasoft.R` reproduces the cross-ancestry meta-analysis summary
(`xGeneMetasoftSummary.pdf`). It is **not** part of Figure 1 and is not run by
`run_all.R`. It requires the MetaSoft Java tool (`METASOFT_JAR`); see the script
header.

## Original source

Reconciled from `QTLAnalysisUtils/notebooks/QTL_paper_Fig1.ipynb` (the authoritative
working notebook). The pre-refactor `Fig1.R` was a linearized export of that
notebook and contained several bugs (undefined assembly variables, a `ptions()`
typo, "casual"→"causal" misspellings, and duplicated read logic), all fixed here.
