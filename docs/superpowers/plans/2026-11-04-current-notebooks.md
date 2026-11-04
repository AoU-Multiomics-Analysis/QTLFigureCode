# Current QTL Notebooks Implementation Plan

> **For agentic workers:** Use superpowers:executing-plans to implement this plan task by task.

**Goal:** Update the five main notebooks and add locally supported supplement notebooks.

**Architecture:** Extend the existing generator with a current-figure index and R scripts.
Keep the original data pipelines. Read local inputs through two explicit root paths.

**Tech Stack:** Python standard library, R, tidyverse, ggplot2, Jupyter notebook JSON.

**Spec:** `docs/superpowers/specs/2026-11-04-current-notebooks-design.md`.

## Global Constraints

- Use the November 2, 2026 submission figures.
- Do not commit controlled inputs or participant rows.
- Preserve the original checkout and source data.
- Keep plot titles and subtitles absent. Use tidyverse R code.
- Distinguish saved-plot rendering from upstream analysis.
- Use a review branch. Do not merge main.

## Review Focus

- Roots with spaces must resolve correctly.
- Missing inputs must produce clear errors before plot code runs.
- Each notebook must regenerate without removing its executed outputs by default.
- S11-S14 must use the EAS revisions and current population order.
- Current figure numbers must replace obsolete main-panel and supplement labels.

### Task 1: Notebook input and generation support

**Files:** `notebooks/current/00_setup.R`, `notebooks/current/index.json`,
`tools/build_current_notebooks.py`, `tests/test_current_notebooks.py`.

**Interfaces:** R setup supplies `paper_file()`, `source_file()`, `output_file()`,
`show_reference()` and `save_plot()`. The index defines each notebook and script.

- [x] Write tests for current notebook paths, explicit roots, missing inputs,
  EAS versions, cell IDs, and preserved execution outputs.
- [x] Run the tests and confirm the initial failure.
- [x] Add setup, index, and generator. Extend the existing generator entry point.
- [x] Run Python tests and parse each R cell.

### Task 2: Main figures

**Files:** `notebooks/current/Fig1.R` through `Fig5.R`; five existing main notebooks.

**Interfaces:** Each script uses Task 1 setup and returns the current plots.

- [x] Add current source and plot steps with accurate method notes.
- [x] Generate and execute the main notebooks from local inputs.
- [x] Check figure values and exact submission reference images.

### Task 3: Supplements and coverage

**Files:** `notebooks/current/supplements/`, `notebooks/supplements/`,
`notebooks/current/coverage.tsv`.

**Interfaces:** Same setup and index as the main figures. One notebook per supported section.

- [x] Add S4, S11-S16, S20-S25, S28, S30 and S31 from verified local sources.
- [x] Record missing, partial, and review sections without fabricating inputs.
- [x] Generate and execute the supplement notebooks. Check counts and output types.

### Task 4: Verification and GitHub delivery

**Files:** `README.md`, current notebook guide, validation report, scoped CI check.

- [x] Run notebook checks, R syntax checks, and suitable repository tests.
- [x] Review all outputs for participant records and personal paths.
- [x] Compare source hashes before and after execution.
- [x] Obtain a separate code review, then fix material findings.
- [ ] Commit named files, push a review branch, and open a pull request.

## Verification record

All 21 notebooks executed: 110 code cells and 45 plot images. The five exact
submission images match their recorded SHA-256 hashes. Input file hashes are
unchanged. Regeneration retained all executed outputs and input records.
Six generator tests, R syntax checks, and the repository R suite passed.
A separate review found two defects in provenance retention and Unicode plotting;
both were fixed and verified. The source-control and clipping caveats are stated
in the notebook guide and section notes. GitHub publication follows this check.
