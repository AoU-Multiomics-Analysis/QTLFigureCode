# Figure 3 cleanup — design

**Date:** 2026-07-23
**Status:** approved
**Branch:** `fig3-cleanup` (off `main`, which now contains the merged Fig1 cleanup, PR #1)
**Supersedes:** `Fig3/Fig3.R` (858 lines, flat notebook export)

---

## 1. Scope

Restructure `Fig3/` into the Fig2-style modular pipeline already applied to Fig1,
extend the manifest-driven input system to a second figure, and fix the defects
found in the notebook export.

Figure 3 covers **cross-QTL colocalization**: FastENLOC colocalization between
each pair of molecular assays (eQTL–pQTL, sQTL–pQTL, eQTL–sQTL) at gene, credible-set
(signal), and variant level, plus the relationship between allele frequency and
cross-assay effect-size correlation.

Unlike Fig1, **nothing is folded in from `QTLAnalysisUtils`** — that package contains
no colocalization code. Every helper Fig3 needs (`LoadEnloc*`, `add_fdr_threshold_label`)
is already inline in `Fig3.R` and stays in `Fig3/00_setup.R`.

---

## 2. Locked decisions

| Decision | Choice |
|---|---|
| Branch base | Branch off `main` after PR #1 merged |
| Loader machinery | **Extract** to `common/manifest.R`, shared by all figures |
| Shared upstream computation | **Memoized derived-data builders** in `00_setup.R`, cached to `fig3/data/` |
| Blank layout slot (panel E) | Placeholder for an analysis **not yet decided**; stub the slot, record in `NOTES.md` |
| `coloc_summary` and other unused analyses | Record intent in `NOTES.md`, drop the code |
| Input provenance | Mark unknown origins local-only (`-`), flag TODOs in the README |
| Theme/palette | Per-figure, **not** shared (Fig1 and Fig3 target very different canvas sizes) |

---

## 3. Target structure

```
common/
└── manifest.R                 # input_manifest / manifest_rows / input_files
                               #   / input_path / load_input   (generic, figure-agnostic)
Fig1/
└── 00_setup.R                 # sources common/manifest.R; otherwise unchanged
Fig3/
├── 00_setup.R                 # libs, FIG3_* paths, font/DPI, theme, SZ_* sizes,
│                              #   enloc loaders, add_fdr_threshold_label,
│                              #   build_or_cache + derived-data builders
├── inputs.tsv                 # manifest: key, type, local_subpath, gs_uri
├── panels/
│   ├── panelA.R               # gene-level: % genes by number of colocalized layers
│   ├── panelB.R               # signal-level: % genes by colocalized signals, + coloc-rate inset
│   ├── panelC.R               # MAF x cross-assay effect-size interaction (forest)
│   ├── panelD.R               # eQTL vs pQTL posterior effect size, faceted by MAF bin
│   └── panelE.R               # reserved slot (stub; content undecided)
├── assemble.R                 # composes panels -> Figure3.pdf / Figure3.png
├── run_all.R                  # panels -> assemble
├── supp_coloc_qc.R            # S1: posterior densities + discovery across thresholds
├── supp_prop_gene_coloc.R     # S2: proportion of genes colocalized, by assay pair
├── data/                      # memoized derived tables (.rds)   [git-ignored]
├── objects/                   # per-panel ggplot objects (.rds)   [git-ignored]
├── NOTES.md                   # deferred analyses + reserved panel E
├── README.md
└── .gitignore
```

`Fig3/Fig3.R` is deleted; this spec plus the README record what it did.

---

## 4. Shared loader extraction

`common/manifest.R` contains, verbatim from `Fig1/00_setup.R`, the generic staging layer:
`.manifest_env`, `input_manifest()`, `manifest_rows()`, `input_files()`, `input_path()`,
`load_input()`. It reads the manifest path from the `MANIFEST` variable, which each
figure's `00_setup.R` sets before sourcing it.

In Fig1, that block is replaced with `source("common/manifest.R")`. Fig1's paths, theme,
palette and folded-in analysis functions are untouched; §4.1 lists the only two other
Fig1 edits this requires.

Rationale: the block is byte-identical across figures and is where the `data.table`
non-standard-evaluation bug was found during the Fig1 work (a manifest column named
`key` shadowing the `key` filter argument). With Fig3–Fig5 to come, a single copy means
a single place to fix.

### 4.1 Generalising the metadata attach

Fig1's `load_input()` hardcodes two metadata columns, mapping manifest `group` → `group`
and manifest `assay` → **`Assay`** (capitalised) on the loaded rows. Fig3 needs a
different metadata column (`type`), so the rule is generalised:

> Every manifest column other than `key`, `local_subpath` and `gs_uri` is metadata, and
> is attached to the loaded rows **under its own name** (skipping `NA` values).

This makes the loader figure-agnostic, but the `assay` → `Assay` rename must move out of
the loader and into the manifest header. Two token changes to Fig1, both required for it
to keep behaving identically:

- `Fig1/inputs.tsv` — header column `assay` → `Assay`.
- `Fig1/panels/panelH.R:8` — `annot$assay` → `annot$Assay` (the only place a Fig1 script
  reads that raw manifest column; every other use reads the attached `Assay`).

With those two changes Fig1's loaded frames are identical to today's, so no Fig1 panel
logic changes. Fig3's `type` column then attaches as `type` with no special-casing —
which also removes nine `mutate(type = 'eQTL-pQTL')`-style calls from the Fig3 code,
since the label now arrives from the manifest.

---

## 5. Derived-data builders

Fig3's panels share a heavy upstream. Under Fig1's strict per-panel loading rule the
AoU annotation table — by far the largest input in the repo — would be read and joined
three times per full build, and the three SuSiE tables parsed five times.

`00_setup.R` therefore provides:

```r
build_or_cache <- function(name, expr, refresh = FALSE)
```

which returns `fig3/data/<name>.rds` if it exists and `refresh` is `FALSE`, otherwise
evaluates `expr`, saves it, and returns it. Each builder wraps one `build_or_cache`
call and takes a `refresh = FALSE` argument that it forwards.

| Builder | Builds | Inputs | Consumers |
|---|---|---|---|
| `susie_finemapping()` | e/s/pQTL SuSiE tables with per-CS `cpip` and variant counts; `cs_list`; gene↔UniProt map; `genes_all_three_assays` | `susie_eqtl`, `susie_sqtl`, `susie_pqtl` | A, B, C, D |
| `coloc_gene_level()` | 3 pairs of `.enloc.gene.out`, FDR-labelled on `GLCP` | `enloc_gene` | A, S1, S2 |
| `coloc_cs_level()` | 3 pairs of `.enloc.sig.out`, FDR-labelled on `RCP` | `enloc_cs` | B, S1 |
| `coloc_variant_level()` | 3 pairs of `.enloc.snp.out`, FDR-labelled on `SCP`, joined to the AoU annotation table | `enloc_variant`, `annotations` | C, D, S1 |
| `coloc_cs_components()` | credible sets collapsed into colocalization components via igraph | `coloc_cs_level()`, `susie_finemapping()` | B |
| `coloc_model_input()` | strand-corrected `max_af`, `log_maf`, per-pair SuSiE posterior means/SDs | `coloc_variant_level()`, `susie_finemapping()` | C, D |

Builders compose by calling each other, so the dependency graph is expressed in code
rather than in a run order. Running `panelA.R` alone never touches the annotation file.

The cache is keyed by name only, not by input content. `NOTES.md` and the README state
that `fig3/data/` must be cleared (or `refresh = TRUE` passed) when an upstream input
changes. This matches the existing `objects/` convention.

---

## 6. Panel decomposition

Layout reproduced from `Fig3.R` lines 772–852:

```
counts <- A + B_with_inset             widths = c(1, 1)
top    <- counts | C                   widths = c(2.5, 1)
bottom <- D | E                        widths = c(1.3, 1)
Fig3   <- (top / bottom) + plot_annotation(tag_levels = "A")
             & theme(plot.tag = element_text(size = SZ_TAG))
```

B's inset (`cs_coloc_rate_plot`) is placed with `inset_element` at
`left = 0.48, bottom = 0.60, right = 0.98, top = 0.98, align_to = "panel"`, unchanged,
and is built inside `panelB.R` so the saved `B.rds` is the composed object.

| Panel | Content | Builders |
|---|---|---|
| **A** | Percent of genes by number of colocalized molecular layers (0–3), restricted to genes tested in all three assays | `coloc_gene_level`, `susie_finemapping` |
| **B** | Percent of genes by number of colocalized signals, with an inset of per-assay credible-set colocalization rate | `coloc_cs_components` |
| **C** | Effect of −log(MAF) on cross-assay effect-size correlation — coefficient and 95% CI per assay pair | `coloc_model_input` |
| **D** | eQTL vs pQTL posterior effect size, faceted by MAF bin (Common / Low-Freq / Rare), with per-bin slope and N | `coloc_model_input` |
| **E** | Reserved slot — content undecided | — |

| Supplemental | Content | Output |
|---|---|---|
| **S1** | GLCP/RCP/SCP posterior densities by pair, over discovery counts across posterior thresholds | `FigS_CrossQTLColocQC.pdf` |
| **S2** | Proportion of genes colocalized, by assay pair | `FigS_CrossQTLPropGeneColoc.pdf` |

Panel E is a stub that returns `ggplot() + theme_void()` and saves it as `E.rds`, so
`assemble.R` needs no conditional. `NOTES.md` records that the slot is reserved,
undecided, and occupies the bottom row at `widths = c(1.3, 1)`.

---

## 7. Defects fixed

### 7.1 Data-correctness — these change figure output

1. **Line 462** — `str_extract(sig_1, "^EN[<43;65;42MSG[0-9]+")` contains a pasted
   terminal escape sequence inside the regex. Under `stringr`/ICU — the engine this
   code uses — the pattern does not compile: ICU reads the inner `[0-9]` as a nested
   character class, leaving the outer `[` unclosed, and `str_extract` raises
   `U_REGEX_MISSING_CLOSE_BRACKET`. **The line therefore cannot execute; the original
   script halts there.** (Under base R's engine the same pattern would instead match
   `ENSG001` by accident, which is likely why the corruption looked survivable.)
   Fix: `"^ENSG[0-9]+"`.

   `[<43;65;42M` is an SGR mouse-click report, so this was almost certainly injected
   into the file by a stray terminal paste **after** the analysis was last run —
   making it file corruption rather than a defect the published figure passed through.
2. **Lines 477–481** — `bind_rows()` lists `cleaned_eqtl_pqtl` twice and never binds
   `cleaned_eqtl_sqtl`. Fix: bind each of the three once. **Of the two, this is the
   defect that could actually have reached published output**, since defect 1 would
   have stopped the script outright.
3. **Line 552** — the annotation `"86.2% shared across ≥1 molecular layer"` is a
   hardcoded string, while the value is computed 130 lines earlier at lines 422–424.
   Fix: compute it from `all_proportions` at plot time.
4. **Lines 681–704 vs 726–753** — the three effect-size models are fit twice. The first
   set omits the `!is.na(...)` filter and uses raw `abs()` where the second uses
   `log1p()`; `allele_frequncy_interaction_model_res` is assigned twice and the first
   value is discarded. Fix: keep only the NA-filtered `log1p` models.
5. **Line 291** — `rare = gvs_max_af < .01` is computed on the uncorrected allele
   frequency and then overwritten at line 307 by the strand-corrected `max_af` version.
   Fix: drop the first.
6. **Line 194** — `QTL_union_gene_list` is a three-way *intersection* (genes present in
   eQTL **and** pQTL **and** sQTL), not a union. Fix: rename `genes_all_three_assays`.

Defects 1 and 2 compound: `cleaned_cs_level_colocs` currently contains only sQTL–pQTL
rows, and that frame builds the igraph edge list behind **panel B and its inset**.
Panel B is expected to change materially once fixed. This is called out in the README
as the result to verify against the previous version.

### 7.2 Output and robustness

7. The assembled main figure (lines 851–852) is built and printed but **never saved**.
   Fix: `ggsave()` `Figure3.pdf` and a `Figure3.png` preview.
8. `ggsave('CrossQTLPropGenesColoc.pdf')` (line 440) and `ggsave('Fig3TopHalf.pdf', …)`
   (line 774) omit `plot=` and therefore save whichever plot was last printed — correct
   in a notebook, fragile under `Rscript`. Fix: pass the plot explicitly.
9. `CrossQTLPosteriorDistributions.pdf` is written twice, at lines 340 and 393, with
   different content; the first is silently overwritten. Fix: one save of the complete
   supplemental.
10. `showtext_opts(dpi = 400)` while every `ggsave()` uses the 300 default, so text is
    mis-scaled. Fix: one `DPI` constant used by both, per the Fig1 convention.
11. Five `options(repr.plot.width = …)` calls have no effect outside a notebook.
    Fix: remove; figure dimensions belong in `ggsave()`.
12. Four `system('gsutil cp <figure> gs://…/PaperFigures/…')` uploads. Fix: drop, per
    the Fig1 precedent; the README documents the upload command instead.

### 7.3 Dead code and duplication

13. `LoadEnlocGeneLevel` is never called, and declares CS-level column names
    (`Signal`, `Num_SNP`, `CPIP_*`, `RCP`, `LCP`) that do not match gene-level enloc
    output — the gene-level files are read with plain `fread`. Fix: delete.
14. `LoadEnlocCSLevel` and `LoadEnlocVariantLevel` are identical but for their column
    names. Fix: one `LoadEnlocTable(path, cols)` applying `type.convert(as.is = TRUE)`,
    which also removes roughly ten scattered `as.numeric()` and
    `across(contains('PIP'), ~as.numeric(.))` calls downstream.
15. Both loaders return only via `setnames()`'s invisible return value. Fix: return
    explicitly.
16. `components(g)` is computed twice (lines 501–502). Fix: compute once.
17. `bind_cols(max_af = fixed_afs$max_af)` joins two frames positionally. Fix: compute
    `max_af` with `mutate()` in place.
18. `add_fdr_threshold_label` stores the threshold in `attr(out, "threshold")`, which the
    downstream dplyr verbs silently drop. Fix: `message()` the threshold per pair so it
    is visible in the run log.
19. Unused/no-op: `rename(cs_id = cs_id)`; `basename()` applied to values already bare
    filenames (lines 148–161); duplicate `library(patchwork)` (line 314); the `type`
    column in `edges` (dropped by the following `distinct(from, to)`); the redundant
    `filter(Coloc == TRUE)` inside panel D's `geom_point` (its data is already filtered).
20. Typos: `effetct_size_comparison`, `CrossQTLDiscoveyrAcrossthresholds`,
    `allele_frequncy_interaction_model_res`.
21. Panels are tagged twice: `p1`, `p2` and `MAF_cross_ome_correlation` each set
    `labs(tag = "A"/"B"/"C")` manually, and `assemble` then applies
    `plot_annotation(tag_levels = "A")` over the top — which is also why
    `cs_coloc_rate_plot` and `inset_plot` need `labs(tag = NULL)` to suppress a tag they
    should never have had. Fix: remove every manual `labs(tag = …)`; tagging is
    `assemble.R`'s job. Per-panel `plot.tag` theming likewise collapses into one
    `& theme(plot.tag = element_text(size = SZ_TAG))` in `assemble.R`, matching Fig1.

---

## 8. Deferred, recorded in `NOTES.md`

**Colocalization covariate analysis.** `coloc_summary` (line 512) builds a per-gene table
joining colocalization rate to median TPM, GeneBayes s_het selection bins, and protein
LOD detection rates — evidently working toward whether cross-assay colocalization rate
tracks expression level, selection pressure, or protein detectability. It is never
plotted or written. The code is dropped; `NOTES.md` records the intent, the joins, and
the three inputs it needs.

Consequently `median_tpms`, `LOD_data` and **`genebayes` all lose every consumer**. The
two `left_join(genebayes, …)` calls at lines 232 and 310 attach `post_mean`, `shet_q` and
`shet_bin`, and no retained panel reads any of those columns. Both joins are therefore
removed. This is a pure simplification and also removes a latent row-duplication risk had
`genebayes` ever contained a repeated `ensg`.

Manifest rows for all three inputs are **kept**, so the provenance record survives; the
README marks them as belonging to the deferred analysis.

Also dropped as unused: `cs_count`, `n_eQTLs`, `n_sQTLs`, `n_pQTLs`, `rare_gene_list`,
and the exploratory `glm(Coloc ~ …)` at lines 719–724.

**Panel E.** Reserved, content undecided.

---

## 9. Input manifest

`Fig3/inputs.tsv` — columns `key`, `type`, `local_subpath`, `gs_uri`
(`-` = not applicable). 10 keys, 16 rows.

`type` holds the assay pair (`eQTL-pQTL`, `sQTL-pQTL`, `eQTL-sQTL`) for the nine enloc
rows and `-` elsewhere. Per §4.1 it is attached to the loaded rows automatically, so the
nine `mutate(type = '…')` calls in the original disappear.

| key | rows | staged to | origin |
|---|---|---|---|
| `susie_eqtl` | 1 | `finemapping/coloc/eQTL_susie.COMB.annotated.cleaned.tsv.gz` | `-` — URI needed |
| `susie_sqtl` | 1 | `finemapping/coloc/sQTL_susie.COMB.annotated.cleaned.tsv.gz` | `-` — URI needed |
| `susie_pqtl` | 1 | `finemapping/coloc/pQTL_susie.COMB.annotated.cleaned.tsv.gz` | `-` — URI needed |
| `enloc_gene` | 3 | `enloc/{eQTL_pQTL,sQTL_pQTL,eQTL_sQTL}.enloc.gene.out` | `-` — URI needed |
| `enloc_cs` | 3 | `enloc/{eQTL_pQTL,sQTL_pQTL,eQTL_sQTL}.enloc.sig.out` | `-` — URI needed |
| `enloc_variant` | 3 | `enloc/{eQTL_pQTL,sQTL_pQTL,eQTL_sQTL}.enloc.snp.out` | `-` — URI needed |
| `annotations` | 1 | `annotated/AoU.COMB.AC10.AN95.biallelic.v9.annotations.tsv.bgz` | `gs://fc-secure-…/checkpoints/` |
| `protein_lod` | 1 | `reference/LOD_proportions_per_proteins.tsv` | `gs://fc-secure-…/` |
| `median_tpm` | 1 | `reference/median_tpms.tsv.gz` | `-` — URI needed; deferred analysis |
| `genebayes_shet` | 1 | `reference/41588_2024_1820_MOESM4_ESM.xlsx` (sheet 2) | `-` — URI needed; deferred analysis |

The 12 enloc and SuSiE rows come from hardcoded paths under
`/home/jupyter/multiomics_qtls/edit/FastEnlocInput/` and bare working-directory
filenames in the original. They are recorded as local-only; a loader failure names the
missing file and the key. The README lists which inputs still need a `gs://` URI.

`enloc/` and `finemapping/` are added to the root `.gitignore` (`finemapping/` is
already there from Fig1).

---

## 10. Out of scope

- Regenerating the FastENLOC outputs themselves — they are treated as inputs.
- Fig4 and Fig5.
- Verifying Fig3 numerically. The inputs are *All of Us* controlled-access data and are
  not available outside the Researcher Workbench, so validation here is limited to parse
  checks, sourcing `00_setup.R`, and loader tests against mock inputs — the same
  strategy used for Fig1. The real test is `Rscript fig3/run_all.R` from the project
  root on Terra.
