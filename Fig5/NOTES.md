# Fig5 — notes

Deliberate differences from the pre-cleanup scripts (`PlotFigure5.R`,
`acvrl1_panels.R`, `spr_browser_panel.R`; preserved during the cleanup in
`_original/`, removed in the final commit). "Result-affecting" items change the
plotted/tested values or how a number is derived; the rest are behavior-preserving.

## Result-affecting (verify against the published figure)

- **Preserved ClinGen AD+AR fan-out.** `variant_effect_by_inheritance`
  (`PlotFigure5.R`, the `left_join(gene_assignments, …)` around L127) joins on a
  gene key that is **not unique**: a gene curated under *both* AD and AR yields two
  rows, so that gene's rare eQTL is counted under **both** inheritance classes in
  Panels B and C. This is the original behavior (it threw a many-to-many warning)
  and is **kept deliberately** to match the published figure. Decided 2026-07-27
  (see `docs/superpowers/specs/2026-07-27-fig5-cleanup-design.md` §2). A dedup to
  one row per (variant, gene) would change the Panel B/C counts.

## Behavior-preserving refactors (no change to the figure)

- **`ref`-column data masking.** The SuSiE table has a column named `ref`, which
  shadows `reference$…` under dplyr data masking. The loaders were split into pure
  seams (`derive_clingen`, `annotate_eqtl`); `annotate_eqtl` takes `clingen_genes`
  as a plain character vector rather than `reference$clingen$…`, so nothing inside
  a dplyr verb collides with `ref`. Same output as the original.

- **Panel B statistic (unchanged, clarified).** Panel B plots
  `mean(|posterior_mean|)` — the **SuSiE effect**, not `log2_aFC` — with a 95% CI
  from `1.96 · sd/√n`. The SE is a genuine standard error and *is* drawn (as the
  CI). Kept exactly as the original; noted only because the column choice is easy
  to misread.

- **`rtracklayer` namespaced.** The SPR browser reads bigWigs via
  `rtracklayer::import()` and never attaches the package, so S4Vectors/BiocGenerics
  cannot mask dplyr verbs (`rename`/`slice`/`desc`). Same tracks, no masking bugs.

- **Tagging moved to `assemble.R`.** The original self-tagged inside each builder
  (e.g. `labs(tag="B")`, and a hand-drawn `"G"` text annotation in the SPR header).
  Panels now return **untagged** objects and `assemble.R` applies the tags: the
  composite panels (A, D, F) are wrapped with `wrap_elements()` so each is one tag
  leaf, and an explicit `tag_levels` sequence `c("A","B","C","D","F","E","G")`
  matches patchwork's depth-first traversal of the row-major layout so the visual
  panels read A–G (plain `tag_levels="A"` would swap E/F, since the bottom 2×2 is
  column-major by meaning: ACVRL1 = D/E on the left, SPR = F/G on the right). The
  zoom insets in D/F are `ignore_tag = TRUE`. Tag glyph **positions** may differ by
  a few points from the original hand-placed tags; panel content is identical. The
  rendered figure is 1,137,605 B vs the pre-cleanup 1,136,396 B (~0.1%, from the
  tag-rendering change only). Pinned by `tests/test_fig5_assembly.R`.

## Dropped

- **`PlotFigureDraft.R`** (2,129 lines) — a dead exploratory draft referenced
  nowhere; not carried over.
- `options(repr.plot.*)` and other notebook-only idioms.
