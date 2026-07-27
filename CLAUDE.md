# CLAUDE.md

Guidance for AI sessions working in this repo. Read this first; it is loaded
automatically at the start of a session.

## What this repo is

R code that generates the figures for a multi-omic **molecular-QTL** manuscript
(expression / protein / splicing QTLs across genetic-ancestry groups in the
*All of Us* cohort). One directory per manuscript figure: `Fig1/` … `Fig5/`, plus
`SuppFig1/`. Shared code lives in `common/`; tests in `tests/`; the notebook
generator in `tools/build_notebooks.py`; design specs and plans under
`docs/superpowers/`.

## Non-negotiables

- **Inputs are *All of Us* controlled-access data. Never commit anything derived
  from a participant** — staged input files, cached `.rds`, or rendered figures.
  `.gitignore` covers the staged dirs, `data/`, `objects/`, and figure outputs.
  **Never `git add -A`** — stage named paths only.
- **Run every script from the project root.** Scripts self-reference as lowercase
  `fig<n>/…` even though the directory is `Fig<N>/` (works on case-insensitive
  macOS / the Workbench; would break on case-sensitive Linux — a known caveat).
- **Tests:** `Rscript tests/run_all.R`. They use mock fixtures — no
  controlled-access data is needed or present locally.

## Cleaning up a figure

The main kind of work here is turning a flat notebook-exported script
(`Fig<N>/Fig<N>.R`, `SuppFig<N>/SuppFig<N>.R`) into the modular, manifest-driven,
tested pipeline the cleaned figures use.

**→ Follow [`docs/figure-cleanup-playbook.md`](docs/figure-cleanup-playbook.md).**
It has the full workflow (brainstorm → design spec → TDD plan → per-task review →
whole-branch review → PR), the target directory structure and conventions
(manifest loader, memoized builders with pure test seams, one `DPI` constant,
pt-vs-mm sizes, tagging only in `assemble.R` + a rendered-tag test, the notebook
generator), a **defect checklist** distilled from real bugs the earlier cleanups
found, and the stacked-branch merge logistics.

**All six figures (Fig1–Fig5, SuppFig1) are now cleaned** and follow the
template, so the playbook applies to *new* figures or to reworking an existing
one. Fig1, Fig3, Fig4, Fig5 and SuppFig1 are the worked examples; Fig2
established the original shape and adopted the manifest last. Read the most
recent `docs/superpowers/specs/*-cleanup-design.md` and
`docs/superpowers/plans/*.md` before starting a new one.

Known gaps, if you are looking for work: Fig1 and Fig2 have no tests, and Fig2,
Fig4 and Fig5 have no `gs://` origins in their manifests (their inputs must be
pre-staged by hand). Each figure's `NOTES.md` lists its own deferred items.

## Conventions at a glance

- Each figure: `00_setup.R` (libs, env-overridable paths, font/DPI, theme + `SZ_*`,
  folded helpers, `build_or_cache` + builders) → `panels/panel<X>.R` (one ggplot
  object each, cached to `objects/`) → `assemble.R` (patchwork → `Figure<N>.pdf`;
  tagging happens here only) → `run_all.R`. Plus `inputs.tsv`, `README.md`,
  `NOTES.md`, and `Fig<N>.ipynb`.
- Any change that alters a **result** goes in that figure's `NOTES.md`, with the
  original line reference and what changed.
- Generating `Fig<N>.ipynb` (via `tools/build_notebooks.py`) is a standard part of
  each cleanup; commit it in the same PR. The tracked notebook must round-trip
  through the generator.
