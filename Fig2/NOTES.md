# Figure 2 — deferred work, open questions, and behaviour changes

Fig2 was the **template** the other figure cleanups were modelled on: it had the
`00_setup.R` → `panels/` → `assemble.R` shape before the shared manifest existed.
It was therefore the last figure to adopt the repo-wide conventions, and the
notes below record that catch-up rather than a from-scratch cleanup.

## Behaviour changes

### 1. Inputs moved to the manifest convention (no result change)

Originally `00_setup.R` defined `read_qtl(f) <- fread(file.path(INPUT, f))` and
every call site passed a **literal filename**, with `INPUT` defaulting to
`input_data`. Fig2 was the only figure with no `inputs.tsv`, so it was also the
only one with no recorded provenance for its inputs.

`inputs.tsv` now declares all six inputs and `read_qtl()` takes a **manifest
key**. Two consequences:

- `INPUT` now defaults to `.` (the staging root) rather than `input_data`, with
  `input_data/` folded into each `local_subpath`. Files resolve to the same
  place, but anyone overriding `FIG2_INPUT` must now point it at the *root*
  that contains `input_data/`, not at `input_data/` itself.
- The optional proteomics input is resolved with
  `tryCatch(input_path("susie_pqtl"), error = ...)` rather than `file.exists()`,
  matching Fig5's optional-input pattern. Behaviour is unchanged: the figure is
  built without pQTL when the file is absent, with a message.

No panel data changed — this is a loading/provenance refactor only.

### 2. Per-panel caches are no longer committed

`Fig2/data/*.rds` and `Fig2/objects/*.rds` used to be tracked in git (~40 MB),
so the notebook's panel cells could run without the raw inputs. They are
**participant-derived** — `data/A.rds` alone held 80,965 variant-level rows with
per-ancestry allele frequencies — which the repo's top-level rule forbids
committing.

They are now git-ignored, matching every other figure. The practical effect is
that `01_prep_data.R` must run before the panels; the notebook's prep cell does
this and its data note says so. Nothing about the computation changed.

> The files were removed from tracking but **remain in git history**. Purging
> them from history requires a separate `git filter-repo` / BFG pass and a
> force-push, which has not been done.

### 3. Output renamed `Fig2.pdf` → `Figure2.pdf`

Every other figure writes `Figure<N>.pdf`; Fig2 wrote `Fig2.pdf`. Renamed for
uniformity, with the root `.gitignore` updated to match. Cosmetic only — if you
have links or scripts pointing at `Fig2.pdf`, they need updating.

## Deferred / open

- **No tests.** Fig2 and Fig1 are the two figures with no coverage under
  `tests/`. The prep step is one long script rather than the pure, separately
  testable builders Fig3–Fig5 use, so there is no seam to test against mock
  fixtures without restructuring `01_prep_data.R` first.
- **`DPI` is 300** here (as in Fig1 and SuppFig1); Fig3–Fig5 use 400. Each
  figure's `ggsave()` DPI matches its own `showtext` DPI, so text scales
  correctly in every case — but the figures are not mutually consistent in
  raster resolution.
- **No `gs://` origins.** All six inputs are locally-derived and carry
  `gs_uri = -`, so they must be pre-staged; the pipeline cannot self-stage them
  on the Workbench. Fig4 and Fig5 have the same gap.

---

## Panel F — the central percent-change bin was mislabelled

**This changes a printed label, not a number.** `Fig2Data/MedianKOSummary.tsv`
labels its central bin `(-10,-10)` — a degenerate interval, since no value can
be both greater than and less than -10. It sits between `(-25,-10]` and
`[10,25)` in the ordering, so it is the **-10% to +10%** bin and the upper
bound lost its sign upstream.

Relabelled to `(-10,10)` in `01_prep_data.R` rather than editing the staged
input, which is controlled-access. `lower` still extracts `-10` from the
string, so the bin's position in `reorder(pct_bin, lower)` is unchanged; only
the axis text differs. The bar itself (median 7,991 genes per individual — the
tallest in the panel) was always plotted correctly.

The underlying generator of `MedianKOSummary.tsv` should be fixed at source;
until then this rename is load-bearing and must not be dropped.

## Panel G — legend moved to the bottom-right corner

The manual legend was anchored at hard-coded `x = 1850 / 2000`, roughly the
middle of the axis. It is now placed at 86% / 89.5% of the data maximum, so it
sits in the bottom-right corner, and the right-hand scale expansion widened
from 0.18 to 0.34 so "Duplication-like" clears the panel edge. The anchors are
computed from the data, so the legend stays in the corner if the counts change.

## Panels E and F — tag placement

E and F share a narrow middle column, and their tags landed against both the
panel content above them and panel D to the left. Rather than buying room with
a large margin -- which also pushed E and F apart -- the tags are moved out to
the top-left corner of each plot with `plot.tag.position = c(0, 1)`, and the
margins are kept small: `margin(t = 26, r = 10, b = 6, l = 24)`, applied with
`&` so it reaches both sub-plots.

**`plot.tag.position` resolves relative to the PANEL region, not the outer
margin.** Widening `l` moves the tag right along with the panel, so it cannot
be used to separate the tag from a y-axis title. Verified empirically at
`l = 140`: the tag stayed locked against E's title.

That matters here because panel E's rotated "Outlier Enrichment" title spans
its whole (short) panel height and ends exactly where F's tag begins. The
clearance therefore has to come from the vertical gap: `t = 26` is the smallest
value at which the F tag clears the title -- at `t = 22` the glyphs touch.
Shortening E's y-axis title is the only way to close the gap further.

## Panel F — the grey bands were ~10x too dark (overplotting, not colour)

The two highlight rectangles used `geom_rect()` **without a `data =` argument**,
so they inherited the plot's data (`median_pc`, 13 rows) and redrew the SAME
rectangle once per row. Thirteen stacked copies at `alpha = 0.05` give an
effective opacity of `1 - 0.95^13 = 0.49`, against panel D's 0.05 — which is
why F's bands read far darker than D's despite both nominally using alpha 0.05.

Replaced with `annotate("rect", ...)`, which draws exactly once, and the fill
changed from `"lightgrey"` to `"grey40"` so it matches
`PlotExpectedExpression()`'s bands exactly. Nothing about the data changed.

**The general trap:** a `geom_*` given constant parameters but no `data =`
inherits the plot's data and draws one copy per row. It is invisible for opaque
fills and compounds for transparent ones.

## Panel G — count labels crowded their points

The labels used `hjust = 1.35` / `hjust = -0.5`, which offsets text by a
fraction of its OWN WIDTH. Short labels therefore sat closer to their point
than long ones: "444" crowded its dot while "1,676" cleared it. Compounding
this, `size = 10` points are several hundred data units wide on this axis, so a
width-proportional gap never reliably cleared the dot's radius.

Both sides now use a constant offset of 6% of the x range, applied as data
columns (`UnderLabelX` / `OverLabelX`) with `hjust = 1` / `0`. Every label has
the same clearance regardless of digit count.

**The offsets are baked into the data on purpose.** Referencing the `LabelPad`
script local from inside `aes()` leaves it lazy: the saved `.rds` keeps only
the symbol and a fresh `assemble.R` session dies with
"object 'LabelPad' not found" — the same trap `panelF.R` documents for
`DelBoxStart` / `DupBoxStart`, and one this change hit before being corrected.
