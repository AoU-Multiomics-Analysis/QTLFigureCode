# Public Manifest Privacy Implementation Plan

> **For agentic workers:** Use superpowers:executing-plans to implement these tasks in the current worktree.

**Goal:** Remove GCS addresses from the public files and keep complete originals in a restricted local folder.

**Architecture:** Save and verify the original manifests before edits. Replace GCS values with `-`, preserve all other manifest fields, and remove copied addresses from documentation. Add a GitHub Actions check for GCS addresses.

**Tech Stack:** Python standard library, TSV manifests, existing R input loader, GitHub Actions.

**Spec:** The user requested removal of GCS addresses and selected a local private copy only.

## Global Constraints

- Keep the complete backup outside Git. Do not upload source data or the backup.
- Keep all manifest rows, columns, local paths, metadata, and public HTTPS sources.
- Keep the `gs_uri` column for the existing R loader. Use `-` for removed values.
- Do not change notebook code, saved plots, or recorded execution hashes.
- Do not rewrite Git history. This change applies to current files.

## Review Focus

- The backup must match each original file and have owner-only access.
- Redaction must change only the GCS values in manifests.
- Public HTTPS references must remain available.
- The check must allow synthetic test inputs and generic scheme descriptions.
- Error messages must name the file and line without printing an address.

### Task 1: Save the private originals

- [x] Save all six manifests and nine documents with GCS addresses outside Git.
- [x] Record the source commit, file sizes, and SHA-256 values in a local inventory.
- [x] Verify exact copies and directory/file modes of 0700/0600.

### Task 2: Redact public files and add a check

**Files:** The six `inputs.tsv` manifests; affected Markdown files; `tools/check_public_manifests.py`; `tests/test_public_manifests.py`; `.github/workflows/public-manifests.yml`.

**Interfaces:** `check_file(path: Path, root: Path) -> list[str]` returns safe diagnostic messages. The CLI checks tracked and unignored files and exits with status 1 if it finds a GCS address. It permits the synthetic `bucket` host only in test scripts.

- [x] Write tests for redacted manifests, HTTPS sources, GCS detection, synthetic fixtures, and safe messages. Run them before implementation.
- [x] Implement the check. Run it against the original files and confirm failure.
- [x] Replace GCS manifest values with `-` and copied document addresses with a placeholder. Update active input instructions.
- [x] Add a small Python-only GitHub Actions job for all pull requests and pushes.
- [x] Compare every public manifest with the private original. Verify that only GCS values changed.
- [x] Run the privacy tests, R manifest tests, and saved-notebook validation.
- [x] Obtain an independent review. No actionable findings were reported.

Publication follows these checks: commit the change, open a pull request, and
verify GitHub Actions. No source data or private backup is included.

Verification: 6 privacy tests, 11 notebook-generator tests, 6 preparation-runner
tests, the R manifest tests, and all 21 saved-notebook checks passed. The private
backup matched all 15 originals. All 94 manifest rows and 3 HTTPS sources were
preserved. Only 42 GCS fields changed. Notebook files and hashes were unchanged.
