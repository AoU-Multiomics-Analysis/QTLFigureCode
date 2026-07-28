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

## Regional window is derived from the fine-mapped variants (panels D and F)

**This changes what panels D and F render.**

`plot_finemapping_with_inset()` used to set the regional x-range from
`range()` of **every tested variant** (`00_setup.R`, originally
`regional_x_range <- range(plot_df$position / 1e6)`), which is the full ~1 Mb
cis window — mostly empty axis either side of the signal. It now takes the span
of the variants at `pip >= pip_thresh` plus `regional_flank_bp` (default 5 kb),
via the pure seam `regional_window_bp()`.

The window is then widened to contain the zoom region and clamped to the tested
range, so the inset's callout box cannot fall outside the panel it points into
and the axis never implies data where none was assessed. With no variant over
the threshold it falls back to every tested variant, as before.

Consequences:

- **Panel D (ACVRL1)** previously auto-derived the full tested range; it now
  renders tight around the credible-set variants.
- **Panel F (SPR)** previously pinned `x_limits = c(72.8817, 72.8934)` on top of
  a `between(position / 1e6, 72.881, 72.8925)` pre-filter — two hand-tuned
  windows that had to stay in sync. The `x_limits` override is gone and the
  window is now derived. **The data pre-filter in `make_spr_browser()` is
  retained** — it selects the locus, not the view — so two magic numbers remain
  there.

`x_limits =` is still honoured as a per-panel escape hatch.

## Gene models come from GENCODE when a GTF is staged

`gene_model_track()` no longer reads the committed tables directly. It resolves
genes for the window through `gene_models_for_window()`, which prefers the
optional `gene_gtf` manifest key and falls back to
`Fig5/reference/gene_{spans,exons}.tsv` when no GTF is staged — so the tests,
the notebook, and anyone without a 1.5 GB GENCODE file still build the figure.
The committed tables are regenerable with `tools/make_gene_models.R`.

GTF-backed extraction keeps protein-coding genes only (`gene_types` in
`common/gene_models.R`); a real window otherwise fills with lncRNAs, antisense
transcripts and pseudogenes and the track becomes an unreadable stack.

**Latent bug fixed in passing:** the track selected genes by **coordinate
alone**, never by chromosome. That was safe only because the two published loci
(chr12 ~51.9 Mb, chr2 ~72.9 Mb) do not overlap numerically — with genes pulled
from a whole-genome GTF, every chromosome's genes in that coordinate range would
have been drawn. `gene_model_track()` now takes `chrom`, and panels D and F pass
`"chr12"` / `"chr2"` explicitly.

Genes whose TSS falls outside the (now tighter) window keep their label: it is
clamped into the panel, and the TSS tick + direction arrow are drawn only when
the TSS is genuinely in view. Previously an out-of-range TSS made ggplot drop
the label silently, leaving an unlabelled gene.

## Panels D/F: markers no longer sliced by the panel edge or the inset

Three separate ways a fine-mapping marker could come out cut in half, all fixed
in `plot_finemapping_with_inset()` (`Fig5/00_setup.R`):

- **The panel's left/right border.** `regional_window_bp()` clamps the window to
  the range of tested variants, so there is *always* a variant sitting exactly on
  each edge, and the x scale draws with `expansion(mult = 0)`. The edge variants
  were rendered as half-diamonds in both D and F. The limits are now padded by
  `edge_pad_x = 1.5%` of the window on each side. The padding is applied to
  `regional_x_limits` rather than to the scale expansion **on purpose** — the
  gene-model track is handed the same limits, so the two panels stay aligned.
  This widens the drawn window slightly (D: 51.870–51.918 → 51.869–51.919 Mb),
  which is a change to the *view*, not to any statistic.
- **The zoom inset's bottom edge.** The inset is drawn on an opaque white
  background, so an edge landing mid-marker slices it. F's PIP-0.06 variant sat
  at 0.056 of the panel height with the inset bottom at 0.05 and showed as an
  arc poking out from underneath. `inset_bottom_clearing_points()` now raises the
  requested bottom to the first height that clears every marker beneath the
  inset's x-span; F settles at ~0.13 and its three low-PIP variants are fully
  drawn, while D (whose variants are all PIP ≥ 0.97) keeps the 0.05 default.
  Nudging by a fixed amount does not work — it just re-slices the next point up
  the PIP ladder, which is what `tests/test_fig5_helpers.R` pins down.
- Note this only moves the *edge*. Variants that fall well under the inset were
  and still are hidden by it; that is the point of overlaying the empty low-PIP
  region.

## Fine-map row (D/F) is taller

`assemble.R` row heights are now written as inches of the 16-inch canvas. D/F
went from ~3.2 in to 4.0 in — they carry a titled regional plot, an inset and a
gene track, and were the most cramped band in the figure. The height comes from
the top row (4.56 → 4.35 in) and the mechanism row E/G (7.52 → 6.99 in).

The two gutters are deliberately **not** equal. The A-C / D-F spacer widens
(0.37 → 0.40 in), but the D-F / E-G one narrows (0.32 → 0.26 in): D/F already
carry an x-axis title along their bottom edge and E/G carry variant captions
above theirs, so the *drawn* separation there is much larger than the spacer
alone and an equal spacer reads as a gap. Figure dimensions are unchanged at
18 × 16 in.
