# Figure Data Preparation Implementation Plan

> **For agentic workers:** Use the existing independent figure tasks and a final review.

**Goal:** Publish the scripts that produce the plotted data from source analysis results.

**Architecture:** Per-figure R preparation scripts share validated local roots,
write to an ignored output tree, and are listed in a dependency manifest.
The notebooks link to these steps and can read their outputs explicitly.

**Tech Stack:** R/tidyverse, Python standard library, GitHub Actions.

**Spec:** `docs/superpowers/specs/2026-11-04-figure-data-preparation-design.md`.

## Global Constraints

- Preserve the November 2 scientific results.
- Do not commit source data or generated data tables.
- Read the original paper and figure checkout in place.
- Do not submit cloud analysis jobs or merge main.
- State when a step still requires an unavailable input or source script.

## Review Focus

- The current novelty rules must use source q values and correct UKBB protein targets.
- Counts must use the same units and denominators as the current figures.
- Dependencies must be ordered, and missing inputs must stop the run clearly.
- New output paths must not overwrite source inputs.
- Notebook links and the preparation overlay must identify actual producing scripts.

### Task 1: Shared preparation support

- [x] Write and run failing tests for dependency order, missing scripts, path roots, and output isolation.
- [x] Add `preparation/common.R`, `tools/prepare_figure_data.py`, and the manifest.
- [x] Test paths with spaces, missing inputs, and unchanged source hashes.

### Task 2: Figure preparation code

- [x] Add current Figure 1 novelty, reference matching, gene categories, and count preparation.
- [x] Add or expose current Figures 2-4 filtering, counts, derived tables, and supplementary summaries.
- [x] Add current Figure 5 selection, ClinVar, ancestry, funnel, and effect preparation.
- [x] Compare regenerated values with current saved inputs; record incomplete chains explicitly.

### Task 3: Notebook connection and checks

- [x] Add notebook preparation links and the explicit prepared-data overlay.
- [x] Regenerate and execute supported notebooks without changing their scientific results.
- [x] Add preparation tests and syntax checks to GitHub Actions.
- [x] Review code and stage named paths only.
- Publish to PR 23. The live check state is recorded on GitHub.

Local verification: all 21 notebooks executed with the prepared-input overlay;
all main reference image hashes match; the full R suite and 14 Python tests pass.
The manifest records the Figure 4A input gap and the partial pre-QC/S30 source boundaries.
