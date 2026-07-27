# Figure 1 — deferred analyses & notes

## Resolved: novel pGene comparison (Panel D)

Panel D shows **novel eGene + sGene + pGene** counts relative to external
consortia. The pGene bar was previously absent: the original notebook referenced
`UKBB_pqtl_path` but never assigned it, so the third bar silently disappeared --
a latent bug rather than a deliberate omission.

It is now wired up. The UKBB-PPP reference is the Nature 2023 supplementary
workbook, declared as `ukbb_pqtl` in `fig1/inputs.tsv` and read by
`LoadUKBBpQTL()` from **sheet 10 (`ST9`), skipping 4 header rows**, filtered to
`cis/trans == "cis"`.

**This changed Panel D.** It gained a third bar; the eGene and sGene bars are
unaffected. Counts on the staged data:

| gene class | eGene | sGene | pGene |
|------------|------:|------:|------:|
| Protein-Coding | 1,207 | 2,035 | 1,783 |
| lncRNA         | 8,696 |   415 |     0 |
| non-coding RNA |   246 |     0 |     0 |

`novel_pQTLs` drops lncRNA explicitly: the proteomics assay cannot measure it,
so a non-zero "novel lncRNA pGene" count would be an artefact of the gene-type
join rather than a finding. The two zeros in the pGene column are therefore
expected, not missing data.

Note this is a **different** UKBB reference from the one on Panel F, which
hard-codes `mean_cs = 9.7` as its credible-set comparison point
(`panels/panelF.R`). That constant is unrelated to this workbook and is still
not derived from data.

## Panel H: annotation source, VEP expansion, and a casing bug

Three compounding changes. **Panel H's result changed** — it went from 3
annotations to 12.

### 1. `annot` now reads the top-level SuSiE tables

The key pointed at per-assay `AnnotateSusie` outputs, one of which
(`annotated/pQTL_susie_annotated.tsv`) had no recorded origin anywhere, so the
panel was permanently unbuildable. It now reads the same three tables Fig2 and
Fig3 use:

    COMB.{expression,splicing,proteomics}_SusieMerged.annotated.ancestry_skew.tsv.gz

All three figures therefore share one SuSiE provenance chain. The `gs_uri` on
these rows was cleared: the old bucket paths point at the superseded
`AnnotateSusie` outputs and would re-stage the wrong files on a clean checkout.

### 2. VEP consequences are expanded before the enrichment

`ComputeFinemappingAnnotationEnrichment()` picks annotation columns up via
`where(is.logical)`. These tables carry consequences as ONE comma-separated
string (`"intron_variant,NMD_transcript_variant"`), ~20% of rows holding more
than one term, so the entire coding / splice / UTR / noncoding-RNA half of
`variant_functional_groups` was invisible: the panel rendered with three cCRE
annotations (PLS, pELS, dELS) and no error.

`.expand_vep_consequences()` in `00_setup.R` one-hots each term into its own
logical column, applied after the `pip` filter so it runs on the reduced set. It
builds the matrix directly rather than `separate_rows()` + `pivot_wider()` —
these tables have ~80 columns and pivot_wider would treat every one as an id
column. It drops VEP's `.` placeholder and refuses to overwrite an existing
column (the cCRE flags are already one-hot).

### 3. Casing bug: both UTR bars were being dropped

`PlotFinemappingEnrichment()` filters against a 12-annotation whitelist by
equality. The whitelist and the tribble spelled the UTR terms
`3_prime_utr_variant` / `5_prime_utr_variant`; VEP emits `3_prime_UTR_variant` /
`5_prime_UTR_variant`.

The two UTR enrichments were computed correctly (3,025 and 1,696 variants,
both significant) and then **silently discarded at the join** — which is exactly
why only 10 of the 12 whitelisted annotations appeared. Names now match VEP.

The whitelist in `annotation_labels` is exactly these twelve:
`splice_donor_variant`, `splice_acceptor_variant`, `splice_region_variant`,
`stop_gained`, `missense_variant`, `inframe_insertion`, `inframe_deletion`,
`5_prime_UTR_variant`, `3_prime_UTR_variant`, `PLS`, `pELS`, `dELS`.

Note there is **no microRNA (or any other noncoding-RNA) entry**, so the
"Noncoding RNA" group defined in `variant_functional_groups`
(`non_coding_transcript_exon_variant`, `non_coding_transcript_variant`,
`mature_miRNA_variant`) is unreachable: those rows are dropped by the
`name %in% names(annotation_labels)` filter and the group can never render.
Widening the panel to noncoding RNA means adding labels to `annotation_labels`,
not just to the tribble.

This bug most likely predates the change: the old `AnnotateSusie` files probably
carried lowercase one-hot columns, so the whitelist matched. Repointing to the
`ancestry_skew` exports, which preserve VEP's native spelling, exposed it.

Both failures share a shape worth remembering: **an annotation whose name does
not match is dropped, never reported.** The panel built cleanly at 3 annotations
and again at 10.

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
