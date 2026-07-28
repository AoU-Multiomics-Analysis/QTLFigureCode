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

- **Inputs are *All of Us* controlled-access data. Never commit individual-level
  data** — staged input files or anything carrying per-participant records (e.g.
  `echo_v4_r2.ancestry_preds.tsv`, keyed by `research_id`). `.gitignore` covers
  the staged dirs, `data/`, `objects/`, and the PDF/PNG figure outputs.
  **Never `git add -A`** — stage named paths only.
- **Executed notebooks are the exception, and it is deliberate.** `Fig<N>.ipynb`
  is committed *with* its rendered panels, so the figures are viewable on
  GitHub without staging 2 GB of inputs. This is judged acceptable because the
  repo is **private**, the rendered panels are aggregate summaries (binned
  counts, per-variant statistics, per-group totals — never per-participant
  rows), and no code cell prints a data frame: every cell ends in a ggplot
  object or a status `message()`. **If you add a cell, keep it that way** — a
  bare data-frame expression at the end of a cell would embed real rows in the
  committed JSON. Re-check before committing an executed notebook.
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

Known gaps, if you are looking for work: Fig2 has no tests and Fig1 has only two
(`tests/test_fig1_enrichment_panel.R` for the enrichment panel and the assay
palette, `tests/test_fig1_enrichment_glm.R` for the standalone PIP-regression
enrichment), and Fig2,
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
  each cleanup; commit it in the same PR. A notebook's **code cells** must
  round-trip through the generator, but an executed notebook no longer matches
  it byte-for-byte, because the generator emits output-free cells. So: change
  the generator → regenerate → re-execute against the staged inputs → commit.
  Regenerating alone silently discards every rendered panel.
