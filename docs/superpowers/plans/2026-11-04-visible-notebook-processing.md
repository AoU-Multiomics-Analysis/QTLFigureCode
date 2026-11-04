# Visible Notebook Processing Implementation Plan

> **For agentic workers:** Use independent figure tasks and one final review.

**Goal:** Show source-result processing and plot construction in all 21 current notebooks.

**Architecture:** An include-aware generator expands existing R transformation
functions and per-notebook processing steps into cells. A local-input runtime
records hashes and keeps intermediate results in memory. Plot builders use
prepared tables in memory and do not read RDS objects.

**Tech Stack:** R/tidyverse, Python standard library, GitHub Actions.

**Spec:** `docs/superpowers/specs/2026-11-04-visible-notebook-processing-design.md`.

## Global Constraints

- Preserve the November 2 scientific selections, counts, models and panel order.
- Do not commit source or generated data tables, or participant-level output.
- Do not read RDS files in generated notebook code or require saved plot templates.
- Show filters, joins, counts and plot construction in visible cells.
- Keep missing-input boundaries explicit. Do not merge or submit cloud jobs.

## Review Focus

- No helper or nested include can hide a saved-object dependency.
- Changes to included code must invalidate saved downstream notebook outputs.
- Fresh notebook sessions must compute from source results; no disk cache shortcut.
- Missing Figure 4A input must produce an explicit skip, not a saved-object fallback.
- Input hashes must include the original source roots and explicit file overrides.

### Task 1: Generator and local runtime

- [x] Add failing include, output invalidation, no-RDS, and memory-only runtime checks.
- [x] Add safe recursive includes and root/file-override tracking.
- [x] Keep notebook code round trips and output privacy checks.

### Task 2: Figure processing and plots

- [x] Rework Figure 2 and S16 from raw result tables with explicit plot builders.
- [x] Rework Figure 3 and S20-S23 from source colocalization and fine-mapping tables.
- [x] Rework Figure 4 and S24/S25 from source TWAS and colocalization tables.
- [x] Rework Figure 1 and S4/S11-S15 from source/reference tables.
- [x] Rework Figure 5 and S28/S30/S31 from source selections and locus inputs.
- [x] Compare regenerated data and model values with accepted inputs locally.

### Task 3: Integration and publication

- [x] Generate and execute all supported notebooks with no RDS inputs.
- [x] Run the full R suite, Python tests, syntax and privacy checks.
- [x] Review the branch, record verification and update the guide.
- Publish a new PR. GitHub records the final check state.

Verification: all 21 notebooks executed from source inputs; 172 code cells and
48 images. No RDS input or prepared overlay was recorded. All source hashes and
main reference-image hashes match. The full R suite, 17 Python tests and eight
inline R fixture files passed. The independent review finding was fixed.
