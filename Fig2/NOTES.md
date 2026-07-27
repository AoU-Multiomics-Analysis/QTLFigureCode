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
