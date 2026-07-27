# Figure-cleanup playbook

How a flat, notebook-exported figure script (`FigN/FigN.R`, `SuppFigN/SuppFigN.R`)
becomes a modular, manifest-driven, tested pipeline in this repo. Written for a
future contributor — human or AI — cleaning up the next figure (Fig5 is the
remaining one). It records both the **workflow** and the **recurring bugs** the
earlier cleanups (Fig1, Fig3, Fig4, SuppFig1) turned up, so you don't rediscover
them.

If you are running with the Superpowers skills, this is a companion to them, not
a replacement: use `brainstorming` → `writing-plans` →
`subagent-driven-development` (or `executing-plans`), and let this doc supply the
repo-specific target shape and the defect checklist.

---

## The workflow

Each figure cleanup is one branch, one spec, one plan, one PR.

1. **Read the whole original script first.** All of it. The flat scripts
   interleave analysis, plotting, and dead/exploratory code; you cannot scope the
   work without reading end to end.
2. **Brainstorm → design spec.** Decide structure, which computations are shared
   (→ builders), which inputs exist, and which defects to fix vs. preserve. Get
   the human to approve. Save to
   `docs/superpowers/specs/YYYY-MM-DD-<figure>-cleanup-design.md` and commit.
   - Ask about anything that changes a *result* and that you cannot infer:
     which of two duplicated computations is authoritative, whether an off-by-one
     was intentional, whether an unused analysis should be kept.
3. **Write the implementation plan.** Bite-sized TDD tasks, complete code in each,
   verbatim from the original where faithful. Save to
   `docs/superpowers/plans/YYYY-MM-DD-<figure>-cleanup.md` and commit.
   - **Self-review the plan against the original.** Every plan so far shipped at
     least one real bug that the self-review or the first task-review caught
     (a broken fixture, an rlang injection that can't work, a fixture whose join
     silently matched nothing). Budget for it.
4. **Execute task-by-task, with a review after each.** A fresh implementer per
   task, then a reviewer that checks spec-compliance *and* fidelity to the
   original. Fix Critical/Important findings and re-review before moving on.
5. **Whole-branch review at the end**, on the most capable model. It catches what
   per-task review cannot — cross-file incoherence, an undocumented result change,
   a stale generated artifact. Every cleanup found something here.
6. **Push and open the PR.** Note results that changed (see `NOTES.md`) and the
   merge order if the branch is stacked (below).

---

## Target structure

Model a new figure on **Fig3** or **Fig4** (the most complete). Supplementals get
the same treatment (see **SuppFig1**).

```
FigN/
├── 00_setup.R        # libraries, env-overridable paths, font/DPI, theme + SZ_*,
│                     #   folded helper functions, build_or_cache + builders
├── inputs.tsv        # manifest: key, [metadata cols], local_subpath, gs_uri
├── panels/
│   └── panelA.R …    # each: source setup, load only its inputs, build ONE
│                     #   ggplot object, saveRDS to objects/, message("Panel X built.")
├── assemble.R        # compose panels (patchwork) -> FigureN.pdf + .png. Tagging lives HERE.
├── run_all.R         # panels -> assemble
├── extra/            # (optional) standalone analyses retained but not in the main figure
├── data/             # memoized derived tables [.gitignored]
├── objects/          # per-panel ggplot objects [.gitignored]
├── FigN.ipynb        # generated notebook (see below)
├── NOTES.md          # deferred work + every deliberate difference from the original
├── README.md
└── .gitignore        # objects/  data/  .DS_Store
```

Run everything from the **project root**. Preserve the original untracked script
in a first commit, then `git rm` it in the last task so it stays recoverable.

---

## Conventions

**Manifest & inputs (`common/manifest.R`).** Each input is a `key` in
`FigN/inputs.tsv` mapping to one or more files: a `local_subpath` (relative to a
staging root that defaults to the cwd) and, where known, a `gs_uri`. The loader
stages a missing file **only if absent**, dispatching on the URI scheme:
`gs://` → `gsutil cp`, `http(s)://` / `file://` → `download.file`. Reach inputs
with `input_path(key)` (single file), `load_input(key)` (row-binds all files for
a key and attaches metadata columns), `input_files(key)` (rows + resolved paths).
Reserved columns are exactly `key`, `local_subpath`, `gs_uri`; any other column is
metadata. Mark unknown origins `-` and list them in the README as needing a URI.

**Memoized builders.** When panels share heavy upstream work, put it behind
`build_or_cache(name, expr, refresh = FALSE)` — it caches to `data/<name>.rds`,
keyed **by name only** (clear `data/` after changing an input; `refresh = TRUE`
does **not** cascade to dependency caches). Panels call only the builders they
need. Keep a **pure seam** inside each builder (e.g. `join_mean_tpms`,
`coloc_cs_edges`, `derive_gene_sets`) that takes frames and returns frames, so the
logic is unit-testable without the real GTF/parquet/network reads.

**Rendering.** One `DPI` constant feeds `showtext_opts()` and every `ggsave()`.
`element_text()` sizes are **points**; `geom_*`/`annotate()`/`geom_label()` sizes
are **millimetres** — never reuse a constant across the two. Font is Inter via
`showtext`, with a `sans` fallback when Google Fonts is unreachable. For a blank
placeholder use `ggplot() + theme_void(base_family = FONT)` (plain `theme_void()`
drops the font and mis-renders any tag on it).

**Tagging happens only in `assemble.R`**, via `plot_annotation(tag_levels = "A")`
or explicit `labs(tag = …)` — never in a panel script. Add a **rendered-tag test**
(`tests/test_<fig>_assembly.R`) that builds the figure and inspects the actual tag
grobs; reading the code is not enough (see the tag bug below). Reuse the
`extract_tag_labels` grid-grob walker from `tests/test_fig3_assembly.R`.

**Notebook (standard deliverable).** Add the figure's config to
`tools/build_notebooks.py`'s `FIGURES` dict (+ `PANEL_NOTES`), then
`python3 tools/build_notebooks.py FigN` and commit `FigN/FigN.ipynb` in the same
PR. The generator pulls code cells from the real pipeline scripts, so notebooks
stay faithful; its setup cell is case-robust (`^(FigN|fign)$`) so it runs on
Terra/Linux. The tracked notebook must round-trip through the generator — after
any generator change, regenerate and commit.

**Tests.** No controlled-access data locally, so: a `make_<fig>_fixture()` in
`tests/fixtures.R` builds mock inputs; the pure seams are tested on those; the
rendered-tag test pins the assembly. Everything else is parse-checked. Install
missing R packages as needed (`igraph`, `MatchIt`, `arrow`, `enrichR`, `hexbin`
have all been needed). Run the whole suite with `Rscript tests/run_all.R`.

**Data policy.** Inputs derive from *All of Us* **controlled-access** data and
must never be committed. `.gitignore` covers the staged input dirs, `data/`,
`objects/`, and the figure outputs. Never `git add -A` — untracked
controlled-access files may sit in the tree; stage named paths only.

---

## Defect checklist — what to hunt for

Every flat script has some of these. The parenthetical says where it bit.

- **Duplicated / dead computations.** A frame computed twice with different
  logic, the second unused; a read whose result is never consumed. *(Fig4:
  `ColocAdjustedConditionalTWAS` computed twice, the first referencing an
  undefined variable; `GeneBayes` read but never used — dropping it removed a
  whole dependency.)* Resolve which is authoritative with the human.
- **Double transforms.** A column transformed (e.g. `log2(x+1)`) and then
  transformed again downstream, while the axis label claims one transform.
  *(SuppFig1 Panel C logged twice.)* Change the plotted/tested values; document it.
- **Hard-coded values that should be computed.** A percentage or count pasted as
  a string while the same value is computed elsewhere. *(Fig3 Panel A: a literal
  "86.2%".)* Compute it.
- **Tag miscounting.** `plot_annotation(tag_levels="A")` letters every *leaf*
  plot. Insets and 2-plot sub-panels each count. *(Fig3: an inset stole a letter
  and shifted C–E to D–F, visibly; SuppFig1: A and D are 2-plot pairs, so the
  figure legitimately has six tags A–F.)* Always verify with a rendered-tag test,
  and use `inset_element(..., ignore_tag = TRUE)` / `wrap_elements()` where
  needed.
- **Unguarded joins that fan out.** A `left_join` on a key that isn't unique
  silently multiplies rows and inflates downstream counts. *(Fig4: annotation and
  posterior joins.)* Assert cardinality (`relationship = "many-to-one"`) with an
  error that names the offending table.
- **Floating-point boundary / NA handling.** Reconstructing a threshold as
  `1-(1-p)` is inexact and can drop the boundary row from its own set; `NA`
  propagating into a label poisons downstream `sum()`/`mean()`. *(Fig3 FDR
  labelling.)* Select by rank on observed values; force missing → `FALSE`.
- **Silent gene-universe / filter changes.** Deriving "all AoU genes" from a table
  that was inner-joined to another cohort narrows the set. *(SuppFig1 detection
  extra.)* Match the original's data source exactly, or document the change.
- **Off-by-one / count changes.** `filter(row_number() < 10)` keeps 9, not 10.
  *(SuppFig1 GO panel.)* Preserve the original count unless the human okays a change.
- **Deprecated / notebook-only idioms.** `stat(x)` → `after_stat(x)`; remove
  `options(repr.plot.*)`; drop unused libraries; add ones actually used but not
  imported (`arrow`, `broom`, …).
- **Ignored parameters, mislabeled stats.** A function that hard-codes a value and
  ignores its own argument; an `SE = sd/n` that is not a standard error and is
  never plotted. Fix or drop, and note it.

Anything you change that affects a **result** goes in `NOTES.md` with the original
line reference and what changed — that file is the record for whoever compares the
new figure to the published one.

---

## Merge & PR logistics

- **The known portability caveat:** scripts self-reference as lowercase `fign/…`
  while the directory is `FigN/`. This works on case-insensitive macOS and the
  Workbench but **fails on case-sensitive Linux**. Consistent across all figures;
  fix repo-wide in one pass if it ever matters.
- **Stacked branches share two files** — `common/manifest.R` and
  `tools/build_notebooks.py`. When several cleanup branches are open, whichever
  merges second hits an add/add conflict on the generator (and a `.gitignore`
  union). Resolution is always **take the later branch's superset generator**
  (it contains the earlier figures' configs plus its own), dedupe `.gitignore`,
  regenerate the notebooks so they round-trip, and confirm `main` is an ancestor
  (`git merge-base --is-ancestor origin/main <branch>`) and
  `git merge-tree --write-tree origin/main <branch>` reports no `CONFLICT`.
- The design spec and plan for each figure live under `docs/superpowers/`; read
  the most recent one as a worked example before starting a new figure.

---

## Status

| Figure | State |
|--------|-------|
| Fig1 | cleaned, merged |
| Fig2 | was already the modular template |
| Fig3 | cleaned, merged |
| Fig4 | cleaned (PR) |
| SuppFig1 | cleaned (PR) |
| Fig5 | empty — next candidate |

Notebooks (`FigN.ipynb`) + the generator are shared; each cleaned figure commits
its own notebook in its PR.
