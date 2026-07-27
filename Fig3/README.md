# Figure 3 — cross-QTL colocalization across molecular layers

Code to generate **Figure 3** (and two supplemental figures) for the manuscript.
The figure asks how often a genetic signal is **shared between molecular layers**
— expression (eQTL), protein (pQTL) and splicing (sQTL) — using FastENLOC
colocalization between each pair of assays at gene, credible-set (signal) and
variant level, and whether that sharing depends on allele frequency.

> **Data note.** The inputs derive from *All of Us* **controlled-access**
> genomic data and are **not** in this repository. See
> [Data & inputs](#data--inputs). Never commit staged inputs or figure
> outputs.

---

## The figure

**Main figure** (`Figure3.pdf` / `Figure3.png`), 18 × 12″:

| Panel | Content | Source |
|-------|---------|--------|
| **A** | Percent of genes by number of colocalized molecular layers (0–3), among genes fine-mapped in all three assays | gene-level FastENLOC + SuSiE |
| **B** | Percent of genes by number of colocalized signals, inset with per-assay credible-set colocalization rate | signal-level FastENLOC + SuSiE |
| **C** | Effect of −log(MAF) on cross-assay effect-size correlation, per assay pair | variant-level FastENLOC + SuSiE |
| **D** | eQTL vs pQTL posterior effect size, faceted by MAF bin | variant-level FastENLOC + SuSiE |
| **E** | *Reserved slot — content undecided (see [`NOTES.md`](NOTES.md))* | — |

**Supplemental figures** (standalone; not built by `run_all.R`):

- `FigS_CrossQTLColocQC.pdf` — GLCP/RCP/SCP posterior distributions by assay
  pair, over the number of colocalizations discovered across posterior
  thresholds.
- `FigS_CrossQTLPropGeneColoc.pdf` — proportion of genes colocalized, by pair.

---

## Repository layout

```
.
├── common/manifest.R            # shared input manifest + gsutil staging (all figures)
├── fig3/
│   ├── 00_setup.R               #   libraries, paths, theme, FastENLOC readers,
│   │                            #     FDR labelling, derived-data builders
│   ├── inputs.tsv               #   input manifest: key, type, local_subpath, gs_uri
│   ├── panels/panelA.R … panelE.R  #  each -> fig3/objects/<X>.rds
│   ├── assemble.R               #   composes panels -> Figure3.pdf / Figure3.png
│   ├── run_all.R                #   panels -> assemble (end-to-end)
│   ├── supp_coloc_qc.R          #   supplemental: posteriors + discovery thresholds
│   ├── supp_prop_gene_coloc.R   #   supplemental: proportion of genes colocalized
│   ├── NOTES.md                 #   deferred analyses, open questions, behaviour changes
│   ├── data/                    #   memoized derived tables (.rds)  [git-ignored]
│   └── objects/                 #   per-panel ggplot objects (.rds) [git-ignored]
├── enloc/ finemapping/ annotated/ reference/   # inputs staged here [git-ignored]
└── Figure3.pdf / Figure3.png
```

### Two caches, and how they differ

**Staged inputs** (`common/manifest.R`) — `inputs.tsv` maps each logical input
(a `key`) to one or more files: a `local_subpath` relative to the staging root
(the cwd by default) and, where known, a `gs://` origin. A file is fetched with
`gsutil` **only when the local copy is missing**, so pre-staged files skip the
copy entirely. Panels reach inputs through `input_path(key)`,
`load_input(key)` and `input_files(key)`.

**Derived tables** (`fig3/data/`) — Fig3's panels share a heavy upstream: the
three SuSiE tables. Six builders in `00_setup.R`
memoize their results, so a panel still declares only what it needs while the
expensive work runs once:

| Builder | Used by |
|---|---|
| `susie_finemapping()` | A, B, C, D |
| `coloc_gene_level()` | A, S1, S2 |
| `coloc_cs_level()` | B, S1 |
| `coloc_variant_level()` | C, D, S1 |
| `coloc_cs_components()` | B |
| `coloc_model_input()` | C, D |

> **This cache is keyed by name, not by content.** After changing an input
> file or a builder, clear `fig3/data/` (or pass `refresh = TRUE`), or you
> will get a stale table back. The same applies to the `objects/` panel
> cache. **`refresh = TRUE` does not cascade**: each builder memoizes under
> its own name and calls its dependencies with their own default
> `refresh = FALSE`, so refreshing one builder can still read a stale cache
> for the tables it depends on. To force a full rebuild, clear the whole
> cache directory rather than refreshing a single builder — see
> [`NOTES.md`](NOTES.md).

---

## Running

Run everything from the **project root**.

```bash
# Full pipeline: build every panel (staging inputs as needed) -> assemble
Rscript fig3/run_all.R

# Supplementals (standalone)
Rscript fig3/supp_coloc_qc.R
Rscript fig3/supp_prop_gene_coloc.R
```

Rebuild a single panel after editing it, then recompose:

```bash
Rscript fig3/panels/panelD.R
Rscript fig3/assemble.R
```

Force a full rebuild of the derived tables:

```bash
rm -rf fig3/data && Rscript fig3/run_all.R
```

Upload finished figures (the scripts no longer do this themselves):

```bash
gsutil cp Figure3.pdf gs://fc-secure-.../PaperFigures/
gsutil cp FigS_*.pdf  gs://fc-secure-.../PaperFigures/SupplementalFigures/CrossQTL/
```

### Tests

`tests/` holds five plain `Rscript` checks that need no controlled-access
data — they run against mock inputs built by `tests/fixtures.R`
(`make_fig3_fixture()` / `use_fig3_fixture()`):

```bash
Rscript tests/run_all.R
```

- `test_manifest.R` — `common/manifest.R`: staging, metadata attach, the
  `key`-column NSE guard.
- `test_fig3_readers.R` — the FastENLOC positional readers and the FDR
  threshold labeller.
- `test_fig3_builders.R` — the six derived-data builders against the mock
  fixture, including `build_or_cache()`'s memoization.
- `test_fig3_edges.R` — regression test for the two defects that reduced the
  credible-set colocalization graph to a single assay pair (see
  [Provenance](#provenance)). `coloc_cs_edges()` is a pure function, so this
  test needs no manifest or fixture.
- `test_fig3_assembly.R` — builds the real panels against the fixture,
  assembles the figure, and inspects the actual grobs patchwork draws to
  assert the rendered panel-letter sequence is exactly A–E (and the QC
  supplemental's is exactly A–B). Added because patchwork's auto-tagging
  silently counted panel B's inset as a taggable sibling, shifting every
  later panel's letter down by one — a bug invisible from reading the panel
  source, since it depends on patchwork's internal tag-walking behavior.

---

## Requirements

- **R** ≥ 4.3
- CRAN: `tidyverse`, `data.table`, `patchwork`, `igraph`, `broom`, `scales`,
  `showtext`, `sysfonts`, `paletteer`

```r
install.packages(c("tidyverse","data.table","patchwork","igraph","broom",
                   "scales","showtext","sysfonts","paletteer"))
```

Text uses **Inter**, fetched from Google Fonts at runtime via `showtext`, with
an automatic fallback to the default sans font. Output DPI is matched to the
`showtext` DPI so text scales correctly.

---

## Data & inputs

Every input is declared in **`fig3/inputs.tsv`**, which is also the record of
each file's origin. Columns: `key`, `type` (assay pair, `-` where not
applicable), `local_subpath` (relative to the staging root, overridable via
`FIG3_INPUT`), `gs_uri` (`-` = local-only).

| Manifest key | Files | Used by |
|---|---|---|
| `susie_eqtl` / `susie_sqtl` / `susie_pqtl` | `finemapping/coloc/{e,s,p}QTL_susie.COMB.annotated.cleaned.tsv.gz` | A, B, C, D |
| `enloc_gene` | `enloc/*.enloc.gene.out` (3 pairs) | A, S1, S2 |
| `enloc_cs` | `enloc/*.enloc.sig.out` (3 pairs) | B, S1 |
| `enloc_variant` | `enloc/*.enloc.snp.out` (3 pairs) | C, D, S1 |
| `protein_lod` / `median_tpm` / `genebayes_shet` | `reference/…` | *(deferred — see NOTES.md)* |

**Only `protein_lod` currently carries a `gs://` URI.** The
rest are local-only and must be pre-staged; a missing one fails with a
message naming the file and the key. See [`NOTES.md`](NOTES.md) for the full
list.

The nine FastENLOC outputs are treated as **inputs**, not pipeline products —
this figure does not run FastENLOC.

To add or re-point an input, edit `fig3/inputs.tsv`; no script changes
needed.

Access to the underlying *All of Us* data is via the *All of Us* Researcher
Workbench under its Data Use Agreement; it is not shared through this
repository.

---

## Provenance

This directory replaces a single ~858-line notebook export (`Fig3.R`). The
restructuring and the defects it fixed are recorded in
[`docs/superpowers/specs/2026-07-23-fig3-cleanup-design.md`](../docs/superpowers/specs/2026-07-23-fig3-cleanup-design.md).
The retired original is recoverable at commit **`dab3c47`** —
`git show dab3c47:Fig3/Fig3.R` — for whoever wants to compare against the
pre-cleanup figure.

**Two defects in the original changed results**, but not in the way they
first appeared to. `Fig3.R:462` held a pasted terminal escape sequence inside
a regex (`"^EN[<43;65;42MSG[0-9]+"`); under `stringr`/ICU — the engine this
code uses — that pattern doesn't compile at all, so `str_extract` raises
`U_REGEX_MISSING_CLOSE_BRACKET` and **the original script halts there**
rather than silently emptying anything. (`[<43;65;42M` is a mouse-click
escape sequence, almost certainly pasted into the file by a stray terminal
event after the analysis last ran — file corruption, not something published
output passed through.) The defect that could actually have reached
published results is the duplicated `bind_rows()` argument at `Fig3.R:477`,
which listed one assay pair twice and never bound eQTL-sQTL — so the
credible-set colocalization graph behind **panel B and its inset** was built
from sQTL-pQTL edges alone. Panel B is expected to differ from the
pre-cleanup version; that difference is the fix, not a regression.

`NOTES.md` also records three further behaviour changes found during the
rewrite that can shift specific rows or panel D's numbers at the margins —
see [`NOTES.md`](NOTES.md) for details.
