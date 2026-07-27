# Figure 4 cleanup — design

**Date:** 2026-07-24
**Status:** approved
**Branch:** `fig4-cleanup` (off `main`, which contains the merged Fig1 + Fig3 cleanups
and the shared `common/manifest.R` and `tests/`)
**Supersedes:** `Fig4/Fig4.R` (774 lines, flat notebook export)

---

## 1. Scope

Restructure `Fig4/` into the Fig2/Fig3-style modular pipeline already applied to
Fig1 and Fig3, extend the manifest-driven input system to a fourth figure, add
memoized derived-data builders for the shared upstream, and fix the defects found
in the notebook export.

Figure 4 covers **colocalization of molecular QTLs with GWAS traits** and the
downstream rare-variant / conditional-TWAS analyses: how many genes colocalize
per trait category and at what credible-set rate, the R² that rare variants add
to TWAS, whether colocalized genes carry inflated conditional-TWAS signal versus
matched controls, and whether that signal enriches in clinically important gene
sets.

As with Fig1/Fig3, **nothing is folded in from `QTLAnalysisUtils`** — every helper
Fig4 needs is already inline in `Fig4.R` and stays in `Fig4/00_setup.R`.

---

## 2. Locked decisions

| Decision | Choice |
|---|---|
| Structure | Match the Fig2/Fig3 template |
| Loader machinery | Reuse the shared `common/manifest.R` |
| Shared upstream computation | **Memoized derived-data builders** in `00_setup.R`, cached to `fig4/data/` |
| `ColocAdjustedConditionalTWAS` duplication | **Keep the second (BH-adjusted) version**, delete the first, rebuild `GeneTWASData` from it |
| Input provenance | Mark unknown origins local-only (`-`), flag TODOs in the README |
| Supplementals | None (Fig4 has none) |
| Tagging | Assigned explicitly in `assemble.R`; verified against rendered grobs |

---

## 3. Target structure

```
Fig4/
├── 00_setup.R                 # libs, FIG4_* paths, font/DPI, theme, SZ_* sizes,
│                              #   palettes (GROUP_LABELS/GROUP_COLORS/RATE_COLUMNS),
│                              #   folded helper functions, build_or_cache + builders
├── inputs.tsv                 # manifest: key, local_subpath, gs_uri
├── panels/
│   ├── panelA.R               # coloc-by-trait-category (counts + CS coloc rate)
│   ├── panelB.R               # R² increase from rare variants (histogram)
│   ├── panelC.R               # conditional-TWAS QQ, colocalized vs matched (blood)
│   └── panelD.R               # gene-set enrichment odds-ratio forest
├── assemble.R                 # composes panels -> Figure4.pdf / Figure4.png
├── run_all.R                  # panels -> assemble (end-to-end)
├── data/                      # memoized derived tables (.rds)   [git-ignored]
├── objects/                   # per-panel ggplot objects (.rds)   [git-ignored]
├── NOTES.md                   # deferred analyses + notes
├── README.md
└── .gitignore
```

`Fig4/Fig4.R` is committed unchanged first (it is currently untracked, so this
preserves the original in history), then removed with `git rm` in the final task.

---

## 4. Folded helpers

All of `Fig4.R`'s inline functions move into `Fig4/00_setup.R` unchanged except
for the bug fixes in §7: `read_any`, `require_columns`, `as_logical_flag`,
`complete_category_colors`, `blend_with_white`, `make_fill`, `reorder_within`,
`scale_y_reordered`, `build_plot_data`, `build_empty_plot`, `build_plot`,
`prepare_qq_match_source`, `calculate_qq_data`, `match_colocalized_pairs`,
`plot_qq`, `fit_enrich`.

The palette / label constants (`GROUP_LABELS`, `GROUP_COLORS`, `RATE_COLUMNS`) and
the blood-trait vector also live in `00_setup.R`.

The library block adds `arrow` and `broom` (used but not imported in the original)
and drops `optparse` **and `readxl`** (imported but unused — see §7.2 item 3a:
`GeneBayes` is read but never consumed). The design-system constants follow the
Fig3 convention: `element_text()` sizes in pt, `geom_*`/`annotate()` sizes in mm,
and one `DPI` constant used by both `showtext_opts()` and every `ggsave()`.

---

## 5. Derived-data builders

`00_setup.R` provides `build_or_cache(name, expr, refresh = FALSE)` (identical to
Fig3's: lazy `expr`, cache to `FIG4_DATA/<name>.rds`, keyed by name only) and the
following builders. Panels call only the builders they need.

| Builder | Builds | Inputs | Consumers |
|---|---|---|---|
| `read_conditional_twas()` | `ConditionalTWAS.cleaned.tsv.gz`, read once | `conditional_twas` | C, D |
| `gene_sets()` | `druggable` / `approved` (drug interactions), Mendelian gene ids (ClinGen Definitive/Strong), heme-immune gene ids (keyword filter on ClinGen disease text) | `drug_interactions`, `clingen` | D |
| `coloc_by_trait()` | `build_plot_data()` over the coloc gene + rate tables (lenient union rate, `min_coloc_genes = 1`) | `coloc_gene`, `coloc_rate` | A |
| `qq_match_source()` | `prepare_qq_match_source()` then `match_colocalized_pairs()` (MatchIt nearest-neighbour propensity match on `log_tpm + log_shet + log_n_rare + abs_vcor`, exact on `coloc_trait`) | `read_conditional_twas()` | C |
| `gene_twas_data()` | the clean `ColocAdjustedConditionalTWAS` (§7.1) then `GeneTWASData` (per-gene summary), joined to `gene_sets()` for the Mendelian / heme-immune / druggable / approved flags | `read_conditional_twas()`, `gene_sets()` | D |

`match_colocalized_pairs()` seeds its RNG (`set.seed(20260713)`), so the cached
match is reproducible; the cache is invalidated by clearing `fig4/data/` (or
`refresh = TRUE`), documented in the README as for Fig3.

C and D share only the raw `ConditionalTWAS` read; their downstream cleaning
differs (C: artifact filter + per-(gene, coloc_trait) slice for the QQ; D:
FDR-adjusted per-`TWAS_path` enrichment table), so each keeps its own builder.

---

## 6. Panel decomposition

Layout reproduced from `Fig4.R` lines 730–761:

```
top_row <- (B | C)                        plot_layout(widths  = c(0.8, 1.2))
top_row <- wrap_elements(full = top_row)  # isolate B/C from D's axis-label spacing
right   <- top_row / D                     plot_layout(heights = c(1, 0.75))
Fig4    <- (A | right)                     plot_layout(widths  = c(7, 8))
```

| Panel | Object | Builders | Content |
|---|---|---|---|
| **A** | `coloc_plot` | `coloc_by_trait` | Per trait category (excluding `neuropsychiatric`, `renal_urological`; `count > 5`; top 10 per category): colocalizing-gene counts as intensity-shaded bars beside the credible-set colocalization rate as a lollipop, with a dashed median reference. |
| **B** | `R2_dist` | *(loads `r2_twas`)* | Histogram of the R² increase (%) that rare variants add, computed as the `adj_R2` difference between AF thresholds `0` and `0.01` for genes with `Pval < 1e-5`. |
| **C** | `p_blood_qq` | `qq_match_source` | Common and rare conditional-TWAS QQ plots over nine blood traits, colocalized genes vs. MatchIt-matched non-colocalized ("Null"), with common/rare display caps. |
| **D** | `p` | `gene_twas_data`, `gene_sets` | Odds-ratio forest (logistic regression, `fit_enrich`): enrichment of rare vs. common conditional-TWAS z-score across four gene-set outcomes — Mendelian (any), Mendelian (immune/hematological), druggable, approved drug target. |

### Tagging

Tags are assigned in `assemble.R` only — no panel script sets `labs(tag = ...)`.
Because the layout is irregular and wraps the B/C row in `wrap_elements`,
`plot_annotation(tag_levels = "A")` traversal is unreliable (the failure mode that
mis-lettered Fig3's inset), so tags are set **explicitly and deterministically**:
`A` on `coloc_plot`, `B` on `R2_dist`, `C` on `p_blood_qq`, `D` on `p`, positioned
with the `& theme(plot.tag = …)` block from the original.

`tests/test_fig4_assembly.R` builds the figure against the fixture, extracts the
rendered tag grobs, and asserts the sequence is exactly `A, B, C, D` on the
intended panels — the rendered-output check added for Fig3, not a code-shape check.

---

## 7. Defects fixed

### 7.1 Data-correctness — these can change figure output

1. **`ColocAdjustedConditionalTWAS` is computed twice** (`Fig4.R` lines 590–612 and
   636–665). The first block references `fdr_rare_sig`, which is never defined
   (`t_fdrboth = fdr_rare_cond < .1 & fdr_rare_sig`), and `GeneTWASData` — which
   Panel D depends on — is built from that first block. The second block recomputes
   the table cleanly (BH adjustment; `t_fdrboth = fdr_rare_cond < 0.10 &
   rare_zscore_fdr < 0.10`) but nothing downstream uses it. Fix: keep the **second**
   version as the single `ColocAdjustedConditionalTWAS`, delete the first, and build
   `GeneTWASData` from it. This removes the undefined-variable reference and is
   recorded in `NOTES.md` as a result-affecting change.
2. **Duplicate flag.** `genuine_rare` and `t_nominal` are the identical expression
   (`p_rare_cond < 0.05 & rare_pvalue < 0.05`). Keep one.

### 7.2 Cleanup — no result change

3. Unused / duplicated reads removed: `gene_dat` (line 328, a duplicate read of the
   coloc gene table), the standalone `gene_data`/`coloc_rate` reads (lines 335–336,
   which duplicate what `build_plot_data` reads internally), and `blank_plot`
   (line 367, never used).
   - **3a. `GeneBayes` is read but never used** (lines 91–94). Every `shet` /
     `log_shet` downstream comes from the `shet` column already present in
     `ConditionalTWAS` (used at lines 400, 623–625), not from the GeneBayes xlsx.
     Removing the read cannot change results. This also drops the `genebayes_shet`
     manifest input and the `readxl` dependency from Fig4 entirely.
4. Library block: add `arrow` and `broom` (used via `arrow::read_parquet` and
   `broom::tidy` but not imported); drop `optparse` and `readxl` (imported, unused —
   `readxl` only served the dead `GeneBayes` read).
5. The CLI-vestige error messages in `build_plot_data` (`"--gene …"`,
   `"--rate_type …"`, `"--coloc_rate …"`, `"--min_coloc_genes …"`) are reworded to
   plain messages naming the argument, since this is no longer an `optparse` CLI.
6. The `heme_immune` leftover comment (`# <- adjust col name if different`) is
   resolved — the column is `diseaseFromSource`, confirmed against the read.
7. `options(repr.plot.*)` (notebook-only) removed; figure size lives in `ggsave()`.
8. One `DPI` constant (400) feeds both `showtext_opts()` and every `ggsave()`.

### 7.3 Deferred / documented

9. `enrichment_results` fits **two** predictor sets (`continuous`, `binary`) but the
   panel plots only the continuous `max_rare` / `max_common` terms. Both are kept
   (the fit is cheap); `NOTES.md` records that the binary-predictor enrichment is
   computed but not plotted.
10. `NOTES.md` records the hand-curated heme/immune keyword regex, and that the
    `GeneBayes` s_het table (`41588_2024_1820_MOESM4_ESM.xlsx`, sheet "Supplementary
    Table 1" — the same file Fig1 (deferred) and Fig3 reference) was read but never
    used, so it is dropped from Fig4; `shet` comes from `ConditionalTWAS`.

---

## 8. Input manifest

`Fig4/inputs.tsv` — columns `key`, `local_subpath`, `gs_uri` (`-` = local-only).
6 keys, 6 rows. Readers other than the default `data.table::fread` are applied in
the builders (the manifest records paths only).

| key | local_subpath | reader | origin |
|---|---|---|---|
| `coloc_gene` | `high_level_coloc_outputs/harmonized/harmonized_coloc.gene.tsv.gz` | `read_any` | extracted from `gs://fc-secure-…/high_level_coloc_outputs.tar.gz` |
| `coloc_rate` | `high_level_coloc_outputs/final_summaries/coloc_rate_by_trait.tsv` | `read_any` | same tar.gz |
| `r2_twas` | `R2_TWAS_combined.tsv` | `fread` | `-` — URI needed |
| `conditional_twas` | `ConditionalTWAS.cleaned.tsv.gz` | `read_any` | `-` — URI needed |
| `drug_interactions` | `reference/interactions.tsv` | `fread` | `-` — URI needed |
| `clingen` | `reference/evidence_clingen/*.parquet` | `arrow::read_parquet` | `-` — URI needed |

The two `high_level_coloc_outputs/` files are extracted from a single tar.gz whose
`gs://` origin the README documents as a manual staging step (`gsutil cp … && tar
-xf …`), mirroring Fig1's external-MetaSoft step. A loader failure names the
missing file and the key.

`/high_level_coloc_outputs/`, `/R2_TWAS_combined.tsv` and
`/ConditionalTWAS.cleaned.tsv.gz` are added to the root `.gitignore` (`/reference/`
is already there from Fig1). The `/data/` and `/susie_local/` protections added on
the `figure-notebooks` branch are re-applied here as well, since those
controlled-access scratch dirs are present in the working tree but that branch is
not yet merged — a defensive measure so nothing controlled-access can be staged.

---

## 9. Testing

The controlled-access inputs are unavailable outside the Workbench, so validation
mirrors Fig1/Fig3:

- `tests/fixtures.R` gains a `make_fig4_fixture()` building a mock Fig4 project
  (small TSV stand-ins for each tabular input; the ClinGen parquet reader is not
  exercised locally because `arrow` may be absent — the gene-set derivation is
  tested on an in-memory frame instead). Fig4 no longer uses `readxl`.
- The pure transforms are tested directly on mock frames: `build_plot_data`
  (rate-type validation, duplicate-trait guard, `min_coloc_genes` trimming),
  `prepare_qq_match_source` (artifact filter, per-(gene, coloc_trait) slice),
  `calculate_qq_data` (rank/expected/observed, capping), the gene-set derivation,
  and `fit_enrich` (a glm on a tiny frame).
- `tests/test_fig4_assembly.R` asserts the rendered A/B/C/D tag sequence.
- Every R file parse-checks; `Rscript tests/run_all.R` stays green.

The real end-to-end run is `Rscript fig4/run_all.R` from the project root on the
Workbench.

---

## 10. Out of scope

- Regenerating the coloc tar.gz or the TWAS / R² inputs — treated as inputs.
- Numeric verification against real *All of Us* data (controlled-access, absent
  locally).
- Fig5 (empty) and any notebook for Fig4 (the R-kernel notebook generator can add
  one afterward, as for Fig1–3).
